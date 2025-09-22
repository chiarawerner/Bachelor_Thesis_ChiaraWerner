This code is adapted from existing code from Aradia Wilms, Alina Podschun and uses Code from Urs Braun (2021) based on the paper Brain network dynamics during working memory
are modulated by dopamine and diminished in schizophrenia (Braun et al., 2021).

First, I selected all ROIs (firstlevel_ROI) and reordered all files (region_reordering) to match the order of CATO files. Second, I performed the NCT analysis (HCP_MEG_full_NCT), whch 
calls for the functions optim_fun and  nct_analysis_task. After that, I performed a repeated measure ANOVA with covariates (all_anovas).

For visualization I used R Studio (table_full) to create a table and Python (global_plots_clean) to create violin plots for all frequency bands and conditions.
