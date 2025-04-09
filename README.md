---
title: "Scripts for MBA paper"
output: markdown
date: "2025-04-09"
---

# Simulation study and analysis scripts for paper "A roadmap for systematic identification and analysis of multiple biases in causal inference"

Authors: Wijesuriya R.\*, Hughes R.A., Carlin J.B, Peters R.L ., Koplin J.J , and Moreno-Betancur M.

\*Author responsible for the code

## Description of the R scripts for the simulation study

-   `1.Params.R` : Runs the data generating models on the HealthNuts case study dataset to store the parameters from these models in the working directory in a csv file named "parameters.csv"

-   `2.sim.gen.R`: Function to simulate the data.

    -   Before generating the data using this script run the scripts 3 and 4 to do the necessary checks

-   `3.Simulation check.R`: Code to sanity check the data generation models

-   `4.Checking power.R`: This script can be used to set the (true) effect size for the average causal effect (ACE)

    -   The value for the true effect is set via trial and error. The regression coefficient for the exposure in the true outcome generation model was set to different values to achieve \~80% power for estimating the marginal ACE using the target analysis models at the conventional 0.05 significance levels. Based on the regression coefficient for the exposure in the true outcome generation model was set to -0.43 and -0.05 in the realistic and enhanced scenarios (also see supplementary documentation of the paper for more details). These final values are saved in "BP_sim.csv"

-   `5. Computing RR and Risk difference_truth.R` : Code to estimate the "true" marginal risk ratio and risk difference by generating a large synthetic population

-   `6.Simulation study_conf.R` : This script implements adjustment for confounding bias for all simulation scenarios and saves the point estimates, standard errors and 95% confidence intervals for the simulation replicates

-   `6.Simulation study_MBA.R` : As above adjusting for all biases (confounding, misclassification of exposure, misclassification of outcome, type 1 selection bias and type II selection bias)

-   `6.Simulation study_misclassX.R` : As above adjusting for misclassification (of exposure) bias

-   `6.Simulation study_misclassY.R` : As above adjusting for misclassification (of outcome) bias

-   `6.Simulation study_type1SB.R` : As above adjusting for type I selection bias

-   `6.Simulation study_type2SB.R` : As above adjusting for type II selection bias

-   `7.Simulation study_results.R` : This script collates all results generated from the adjustment approaches and estimates the performance measures

    -   The performance measures (bias, relative bias, empirical standard error, model-based standard error, coverage and bias-eliminated coverage) are estimated using the `rsimsum` package.

    -   The final results are saved as .xlsx files and named as *"effectmeasure.results.simulationscenrio.xlsx"* For example, results for the risk difference under enhanced scenario with correct bias parameters : *RD.results.enhanced1.xlsx ;* results for the risk ratio under realistic scenario with misspecified bias parameters : *RR.results.realistic2.xlsx*

## Description of the R scripts for the case study

-   `Case study application_imputation approach.R` : Applies each of the one-at-a-time adjustments and the simultaneous adjustments to the real case study data and generate Table S6 and Figure 4

    -   Note: Case study data not provided due to to ethics requirements but are available from study investigators upon reasonable request. Please direct any case study data requests to health.nuts\@mcri.edu.au.
