# Inference for Coherent Systems with Weibull Distributed Component Lifetimes
## Authors
Pallavi N. Pawar
(Department of Statistics, S.P. College, Pune)
## Co-author
Dr. Madhuri Kulkarni
(Department of Statistics, Savitribai Phule Pune University)
## Description
This repository contains the data and R code for the research work on statistical inference for coherent systems when component lifetimes follow a Weibull distribution. The study applies maximum likelihood estimation (MLE) to both real data and simulated data to analyze system reliability.

## Data Description
The real dataset consists of failure times of air conditioning systems installed in aircraft. This dataset is used to illustrate the proposed methodology and assess its practical applicability.

## Methodology
- Weibull distribution is used to model component lifetimes  
- Conditional Maximum Likelihood Estimation (MLE) is employed for parameter estimation  
- Inference is carried out for coherent systems
- A simulation study is conducted to validate the performance of the estimators  

## Repository Structure
- `data/` : Aircraft failure time dataset  
- `code/` :  
  - `aircraft_data_analysis.R` : R code for real data analysis  
  - `simulation_study.R` : R code for simulation study
  -All results reported in the paper can be reproduced using the provided R scripts.

## How to Reproduce the Results
1. Open R or RStudio  
2. Run `aircraft_data_analysis.R` for real data results  
3. Run `simulation_study.R` for simulation results  
4. Ensure the dataset in the `data/` folder is properly loaded  

## Requirements
- R (version 4.0 or above recommended)  
- Required R packages should be installed before running the scripts  

## Note
This repository is created to ensure reproducibility of the research work. The provided data and code allow replication of all results presented in the paper.
