# SMART-D data analysis

A set of scripts used to analyze SMART-D clinical trial data with Drs Boyu Ren, Lauren Borchers and Manuel Kuhn.

*smart_analysis_final*: Analysis of the EMBARC and SMART-D data. The first part of the script uses EMBARC data to derive predictive algorithms for treatment response to sertraline and bupropion, which form the basis for the two biomarkers used in the SMART-D study. We considered three candidate models: logistic regression, logistic regression with lasso penalty, and random forest, and selected the one that maximized the cross-validated AUROC. An optimal cut-off point was identified in the same process as the point that maximized the distance from (sensitivity, specificity) to (1,1). The final algorithm was built by fitting the selected model on the full data and applying the optimal cut-off. The second part of the script contains the SMART-D response analysis using linear mixed models, fit with the lmer function in the R package lme4. Four models were fit to compare treatment response: (1) between sertraline and bupropion; (2) across three patient groups defined by sertraline and bupropion biomarker status (double positive, single positive, double negative); (3) between patients whose treatment was consistent with their biomarker status and those whose treatment was not; and (4) across all four biomarker status groups. Sex, age, time, group, the interaction between time and group are included as the fixed effects. A random intercept and a random slope for time are included as the random effects.

*smart_emabarc_nature_mental_health.Rmd* This script analyzes SMART and EMBARC study data. Main analyses examine longitudinal changes in depression (MADRS) across biomarker groups, antidepressant medication (bupropion vs. sertraline), and consistent medication assignment.

Code includes:
• Combining, merging, and cleaning dfs
• Longitudinal analysis (linear mixed-effects models)
• Descriptive statistics and group comparisons (t-tests, chi-squares)
• Figures included in the manuscript
• Most statistics included in the manuscript 
