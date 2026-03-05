rm(list = ls())
set.seed(260225)

# install.packages("survival")
library(survival)



# ============================================================
# 1. Set the True parameters
# ============================================================



# True coefficients for the treatment assignment model
# logit P(T = 1 | age, sex, com1, com2) = X * beta_ps_true (+ site shift if used)
beta_ps_true   <- c(-0.3, 0.04, 0.6, 0.8, -0.5)

# True coefficients for the survival outcome model
# h(t|T, Z ) = h0(t) * exp(beta_cox_true^T * (T, age, sex, com1, com2))
beta_cox_true <- c(-0.2, 0.02, 0.3, 0.5, -0.2)

# baseline hazard
# cox_lambda_0 used only for data generation
cox_lambda_0 <- 0.002

# Administrative censoring time: min(T_event, C_random, c_max)
c_max   <- 5 * 365



# ============================================================
# 2. Data generation per site
# ============================================================



# create specific parameters to make heterogeneity in K sites
# different age distribution(mean, sd)
# different sex and comorbidities(질병 유무) distribution
# create site-specific intercept
make_site_profiles <- function(K) {
  data.frame(
    site = 1:K,
    age_mu = runif(K, 55, 75),
    age_sd = runif(K, 2, 5),
    p_sex1 = runif(K, 0.35, 0.65),
    p_com1 = runif(K, 0.10, 0.45), # Comorbidity1
    p_com2 = runif(K, 0.05, 0.35), # Comorbidity2
    ps_shift = rnorm(K, 0, 0.5)
  )
}



# Sigmoid function
sigmoid <- function(x) 1 / (1 + exp(-x))



generate_one_site <- function(n, prof) {
  
  # Generate baseline covariates within one site
  age_raw <- rnorm(n, prof$age_mu, prof$age_sd)
  sex  <- rbinom(n, 1, prof$p_sex1)
  com1 <- rbinom(n, 1, prof$p_com1)
  com2 <- rbinom(n, 1, prof$p_com2)
  
  # Scaling for analysis stability
  age  <- (age_raw - 65)/10
  
  # Treatment assignment via logistic propensity score model
  # logit P(T = 1 | Z) = Z_ps %*% beta_ps_true + site_shift
  Z_ps <- cbind(1, age, sex, com1, com2)
  ps <- sigmoid(Z_ps %*% beta_ps_true + prof$ps_shift)
  treat <- rbinom(n,1,ps)
  
  # Survival outcome generation under Cox
  # hazard: h(t | T, Z) = cox_lambda_0 * exp(eta), h0(t) = cox_lambda_0
  Z_surv <- cbind(treat, age, sex, com1, com2)
  eta <- Z_surv %*% beta_cox_true
  
  # Inverse transform sampling for exponential baseline
  # t_event = -log(u) / (cox_lambda_0 * exp(eta))
  u <- runif(n)
  t_event <- -log(u)/(cox_lambda_0*exp(eta))

  # Censoring(random censoring, administrative censoring)
  c_rand <- runif(n,0,c_max)
  time  <- round(pmin(t_event,c_rand,c_max))
  event <- as.integer(t_event <= pmin(c_rand,c_max))
  
  data.frame(site = prof$site, id = 1:n,
             age, sex, com1, com2,
             treat, time, event)
}



# Merge data of all site
generate_multisite <- function(site_sizes){
  profs <- make_site_profiles(length(site_sizes))
  sites <- lapply(1:length(site_sizes),
                  function(k) generate_one_site(site_sizes[k], profs[k,]))
  list(sites=sites,
       pooled=do.call(rbind, sites))
}



# Make samples of ten sites
K <- 10
site_sizes <- sample(200:800, K, replace = TRUE)
sim <- generate_multisite(site_sizes)



# ============================================================
# 3. Federated Propensity Score
# ============================================================



# Compute site-specific score vector and Hessian matrix for logistic regression (ps model)
score_hessian_site_ps <- function(data, beta){
  X <- cbind(1, data$age, data$sex, data$com1, data$com2)
  y <- data$treat
  p <- sigmoid(X%*%beta)
  score <- t(X)%*%(y-p)
  W <- as.vector(p*(1-p))
  hess <- -t(X)%*%(X*W)
  list(score = score, hess = hess)
}



# Newton-Rapson Algorithm for federated logistic regression
# Aggregates score and hessian of all sites
update_ps_newton <- function(sites, beta_init){
  beta <- beta_init
  for(i in 1:50){
    U <- 0; H <- 0
    
    # Sum score and hessian over all sites
    for(k in sites){
      sh <- score_hessian_site_ps(k,beta)
      U <- U + sh$score
      H <- H + sh$hess
    }
    
    # Newton-Rapson update
    # beta_new = beta - H^{-1} U
    beta_new <- beta - solve(H,U)
    
    # Convergence criterion
    if(max(abs(beta_new - beta)) < 1e-10) break
    beta <- beta_new
  }
  beta
}



# Estimate PS coefficients using federated Newton-Rapson Algorithm
beta_ps_fed <- update_ps_newton(sim$sites,rep(0,5))



# Logistic regression using pooled data
fit_ps_pooled <- glm(treat ~ age + sex + com1 + com2,
                   data = sim$pooled, family = binomial())
beta_ps_pooled <- coef(fit_ps_pooled)



# ============================================================
# 4. IPTW
# ============================================================



# Fit logistic regression(for Propensity Score)
predict_ps <- function(df,beta){
  X <- cbind(1, df$age, df$sex, df$com1, df$com2)
  as.vector(sigmoid(X %*% beta))
}
sim$pooled$ps_hat <- predict_ps(sim$pooled, beta_ps_fed)
p_treat <- mean(sim$pooled$treat)



# Compute stabilized inverse probability of treatment weights(IPTW)
sim$pooled$sw <- ifelse(sim$pooled$treat == 1,
                        p_treat / sim$pooled$ps_hat,
                        (1 - p_treat) / (1 - sim$pooled$ps_hat))



# Attach weights to sites
w_map <- sim$pooled[, c("site", "id", "sw")]
sim$sites <- lapply(sim$sites, function(df)
  merge(df, w_map, by = c("site", "id"), sort = FALSE))



# ============================================================
# 5. Federated Cox (Sharing score/information)
# ============================================================



# Calculate score and hessian for using Newton-Rapson
site_cox_summary <- function(df, beta, event_times){
  time <- df$time
  event <- df$event
  x <- df$treat
  w <- df$sw
  exp_bx <- exp(beta * x)
  
  z0 <- w * exp_bx
  z1 <- w * x * exp_bx
  z2 <- w * (x^2) * exp_bx
  
  ord <- order(time)
  time_s <- time[ord]
  z0_ord <- z0[ord]
  z1_ord <- z1[ord]
  z2_ord <- z2[ord]
  
  # rearrange data because of time event
  S0_sum_i_to_end <- rev(cumsum(rev(z0_ord)))
  S1_sum_i_to_end <- rev(cumsum(rev(z1_ord)))
  S2_sum_i_to_end <- rev(cumsum(rev(z2_ord)))
  
  m <- length(event_times)
  d <- d1 <- S0 <- S1 <- S2 <- numeric(m)
  
  ev_idx <- which(event == 1)
  if(length(ev_idx) > 0){
    d_by_t <- tapply(w[ev_idx], time[ev_idx], sum)
    d1_by_t <- tapply(w[ev_idx] * x[ev_idx], time[ev_idx], sum)
    key <- as.character(event_times)
    d[!is.na(d_by_t[key])] <- d_by_t[key][!is.na(d_by_t[key])]
    d1[!is.na(d1_by_t[key])] <- d1_by_t[key][!is.na(d1_by_t[key])]
  }
  
  idx <- findInterval(event_times - 0.5, time_s) + 1
  in_range <- idx <= length(time_s)
  S0[in_range] <- S0_sum_i_to_end[idx[in_range]]
  S1[in_range] <- S1_sum_i_to_end[idx[in_range]]
  S2[in_range] <- S2_sum_i_to_end[idx[in_range]]
  
  list(d = d, d1 = d1, S0 = S0, S1 = S1, S2 = S2)
}



# Newton-Rapson Algorithm
federated_cox_newton <- function(sites, pooled){
  event_times <- sort(unique(pooled$time[pooled$event == 1]))
  beta <- 0
  for(it in 1:50){
    d_tot <- d1_tot <- S0_tot <- S1_tot <- S2_tot <- numeric(length(event_times))
    for(df in sites){
      s <- site_cox_summary(df, beta, event_times)
      d_tot <- d_tot + s$d
      d1_tot <- d1_tot + s$d1
      S0_tot <- S0_tot + s$S0
      S1_tot <- S1_tot + s$S1
      S2_tot <- S2_tot + s$S2
    }
    ratio1 <- S1_tot / pmax(S0_tot, 1e-12)
    ratio2 <- S2_tot / pmax(S0_tot, 1e-12)
    U <- sum(d1_tot - d_tot * ratio1)
    I <- sum(d_tot * (ratio2 - ratio1^2))
    beta_new <- beta + U / I
    if(abs(beta_new - beta) < 1e-10) break
    beta <- beta_new
  }
  
  I_all <- I
  se_model <- sqrt(1 / I_all)
  
  Uks<-numeric(length(sites))
  for(k in 1:length(sites)){
    s <- site_cox_summary(sites[[k]], beta,event_times)
    ratio1_k <- s$S1/pmax(s$S0, 1e-12)
    Uks[k] <- sum(s$d1 - s$d * ratio1_k)
  }
  se_rob <- sqrt((1 / I_all)^2 * sum(Uks^2))
  
  list(beta = beta, se_model = se_model, se_rob = se_rob)
}


# Estimate beta using Newton-Rapson algorithm
fed_cox <- federated_cox_newton(sim$sites, sim$pooled)



# ============================================================
# 6. Pooled Cox
# ============================================================



# Fit cox regression(pooled data)
fit_pooled<-coxph(Surv(time, event) ~ treat,
                data = sim$pooled,
                weights = sim$pooled$sw,
                ties = "breslow")

beta_pooled <- coef(fit_pooled)[["treat"]]
se_pooled <- sqrt(vcov(fit_pooled)["treat","treat"])



# ============================================================
# 7. CI comparison and Visualization
# ============================================================



# Calculate CI
out <- data.frame(
  Method = c("Pooled", "Federated"),
  Beta = c(beta_pooled, fed_cox$beta),
  HR = exp(c(beta_pooled, fed_cox$beta)),
  CI_L = exp(c(beta_pooled - 1.96 * se_pooled,
             fed_cox$beta - 1.96 * fed_cox$se_model)),
  CI_U = exp(c(beta_pooled + 1.96 * se_pooled,
             fed_cox$beta + 1.96 * fed_cox$se_model))
)

print(out, row.names = FALSE)

cat("\nAbs diff beta:",
    abs(beta_pooled - fed_cox$beta), "\n")
cat("True HR:",
    exp(beta_cox_true[1]), "\n")



# Visualization
# install.packages("forestplot")
library(forestplot)

table_text <- cbind(
  c("Method", out$Method),
  c("aHR", sprintf("%.2f", out$HR)),
  c("95% CI", sprintf("(%.2f, %.2f)", out$CI_L, out$CI_U))
)

hr_values <- c(NA, out$HR)
ci_lower  <- c(NA, out$CI_L)
ci_upper  <- c(NA, out$CI_U)

forestplot(
  labeltext = table_text,
  mean = hr_values,
  lower = ci_lower,
  upper = ci_upper,
  zero = 1,
  is.summary = c(TRUE, rep(FALSE, nrow(out))),
  boxsize = 0.2,
  lineheight = unit(1, "cm"),
  col = fpColors(box = "black", line = "black", zero = "gray50"),
  lwd.ci = 2,
  vertices = TRUE,
  xlab = "adjusted Hazard Ratio",
  graph.pos = 2
)
