# Inference for Coherent Systems with Weibull Distributed Component Lifetimes
## Authors
Pallavi N. Pawar
(Department of Statistics, S.P. College, Pune)
## Co-author
Dr. Madhuri Kulkarni
(Department of Statistics, Savitribai Phule Pune University)
## Description
This repository contains the dataset and R code used in the research on inference for coherent systems with Weibull-distributed component lifetimes. The study applies conditional maximum likelihood estimation (MLE) methods. 
## Data Description
The dataset consists of failure times of air conditioning systems installed in aircraft. It is used to illustrate the proposed methodology and to assess its practical applicability.
The dataset is taken from:
Proschan, F. (1963). *Theoretical explanation of observed decreasing failure rate*. Technometrics, 5(3), 375–383.

## Methodology
- The failure times of the air conditioning systems are modeled using a 2-out-of-4:F system.
- The Weibull distribution is used to model component lifetimes.
- Conditional Maximum Likelihood Estimation (MLE) is employed for parameter estimation.
- Statistical inference is carried out for coherent systems.
- A simulation study is conducted to evaluate the performance of the proposed estimators.

## Repository Structure
- `data/` : Aircraft failure time dataset  
- `code/` :  
- `aircraft_data_analysis.R` : R code for modelling the failure times of the Air conditioning systems
- `simulation_study.R` : R code for comparison with existing estimators through simulation study
- All results reported in the paper can be reproduced using the provided R scripts.

## How to Reproduce the Results
1. Open R or RStudio  
2. Run `aircraft_data_analysis.R` for real data results  
4. Ensure the dataset in the `data/` folder is properly loaded
5. Run `simulation_study.R`for simulation results

## Requirements
- R (version 4.0 or above recommended)  
- Required R packages should be installed before running the scripts  

## Note
This repository is created to ensure reproducibility of the research work. The provided data and code allow replication of all results presented in the paper.
