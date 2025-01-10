Code for paper "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations"

Authors: Mareike Buhl, Eugen Kludt, Lena Schell-Majoor, Paul Avan, Marta Campi
 
Submitted to International Journal of Audiology (10.01.2025)
Preprint: 

The following Matlab scripts contain the code used in the paper. The scripts r2_s1_preprocessing, r2_s2_sii_calc_opt, r2_s3_sii_combine, and r2_s4_srt_assumptions_v3 were conducted via remote analysis on servers of Hanover Medical School. The remaining scripts were conducted locally, based on the results of the remote analysis. 
Results are contained for scripts s4 and s5, which are needed for the local analysis. 

## 1) Remote: Preprocessing  
r2_s1_preprocessing.m
- Exclude patients without Freiburger results or audiogram 
- Audiogram interpolation 
- Better ear 
- Bisgaard 
- Check for unique audiograms to save computing time in remote analysis (s2) 

## 2) Remote: Calculate SII 
r2_s2_sii_calc_opt.m 
- Load T_AG_SII estimated in s1 (only unique audiograms) 

## 3) Remote: Combine SII results to one table 
r2_s3_sii_combine.m
- For SII results (per level) and slope (per patient)

## 4) Remote: Calculate SRTs by linear fit (different versions) for cases based on data availability 
r2_s4_srt_estimation.m   
- Load T_model_indiv for individual SII slopes (only unique audiograms)
- Join the SII slopes with preprocessed table T_all from s1  

## 5) Local: Additional calculations 
r2l_s5_analysis.m
- Plomp analysis 
- Error calculation 
- Prediction of SRT difference 
- Statistical analysis 

## 6) Local: Plot scripts 
r2l_s6_plots_fig2.m
r2l_s6_plots_fig3.m
r2l_s6_plots_fig4.m
r2l_s6_plots_fig5.m  
r2l_s6_plots_fig6.m
r2l_s6_plots_figA1.m

Plot scripts corresponding to the respective figure number. 







