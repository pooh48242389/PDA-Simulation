About this project \\
I created this simulation with dummy data to follow the paper "Federated target trial emulation using distributed observational data for treatment effect estimation".

What I did
Heterogeneous Data: Since real-world hospitals have different patient data, I generated the dummy data to have different covariate distributions for each institution.
Federated Learning: The institutions only share their score and hessian values to update the models. Raw data is not shared.
IPTW & Models: I used a logistic regression model to get propensity scores, applied IPTW (Inverse Probability of Treatment Weighting), and then fitted a Cox regression model.
Algorithm Comparison: I tested two different algorithms to update the models:

FL-TTE_FedProx_Simulation.R: The method used in the original paper.
FL-TTE_NewtonRapson_Simulation.R: The algorithm I applied for this simulation.

Result
Based on the simulation, the Newton-Raphson algorithm gave hazard estimates that are slightly closer to the pooled data results compared to the FedProx method.
