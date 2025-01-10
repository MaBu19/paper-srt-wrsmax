% Code for paper
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations"
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology
%
% Local script 5: Combination of SII results to one table
%
% Mareike Buhl
% mareike.buhl@pasteur.fr
%

clear all;
close all;
clc;

% change here to have all paths correct:
dataflag = 1; % 1: remote analysis (EK), 2: synthetic example analysis (MB)

warning('off');

%% read data
if dataflag == 1
    addpath(genpath('./functions_remote/')); % all final functions need to be put here
    res_folder = './results_remote/';
    fig_folder = './figures_remote/';
elseif dataflag == 2
    path_str = '_remote'; % '_remote' | ''
    addpath(genpath(['./functions' path_str '/']));
    res_folder = './results_example/';
    fig_folder = './figures_example/';
end

pres_paper = 'paper';

if ~exist([res_folder filesep pres_paper], 'dir')
    mkdir([res_folder filesep pres_paper]);
end

sflag = 1;

%% load data (remote analysis)

load([res_folder pres_paper filesep 's4_T_fit']);
T_fit_emp_all_sii = T_fit; 


% reference data
srt_ref = 29.3; % NH reference FMST in quiet

% SII error
sii_error = 0.00084;

%% Plot properties

% color-code for cases
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case == 4) = 1;
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case == 3) = 2;
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 1) = 3;
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2) = 4;
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case == 5) = 5;
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case == 6) = 6;
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case == 0) = 7;
T_fit_emp_all_sii.case_color(T_fit_emp_all_sii.case <=1) = 10; % dummy

% remove data points that don't belong to any of the
% cases, therefore still have case_color = 0
T_fit_emp_all_sii(T_fit_emp_all_sii.case_color==0,:) = [];
% remove patients without PTA
T_fit_emp_all_sii(isnan(T_fit_emp_all_sii.pta_spl),:) = [];

%% Additional calculations (not required to be performed as remote analysis)

thr_c = -10; % -10; % threshold to define what is still considered as consistent

% estimate non consistent SI and PTA data points
T_fit_emp_all_sii.srt_pta_diff_sii = T_fit_emp_all_sii.fsii_indiv_srt-T_fit_emp_all_sii.pta_spl; % >= 0: SRT should be "audible" / PTA is trusted
T_fit_emp_all_sii.srt_pta_diff_dd = T_fit_emp_all_sii.fall_srt-T_fit_emp_all_sii.pta_spl;
T_fit_emp_all_sii.srt_pta_diff_80 = T_fit_emp_all_sii.f1_srt-T_fit_emp_all_sii.pta_spl;

% report some numbers -> flowchart in paper:
% sum(T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c)
% sum(T_fit_emp_all_sii.srt_pta_diff_sii < thr_c) % together
% sum(~isnan(T_fit_emp_all_sii.fsii_indiv_srt)), i.e. only the patients for
% which an SII-slope was available
% sum(T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c)
% sum(T_fit_emp_all_sii.srt_pta_diff_dd < thr_c)
% sum(T_fit_emp_all_sii.srt_pta_diff_80 >= thr_c)
% sum(T_fit_emp_all_sii.srt_pta_diff_80 < thr_c)

% additionally exclude patients for which estimated SRT (linear part of
% psyfun with assumed SII slope) is not consistent with maxSI (dBopt)
T_fit_emp_all_sii.SI_sii = T_fit_emp_all_sii.fsii_indiv_intercept + T_fit_emp_all_sii.fsii_indiv_slope.*T_fit_emp_all_sii.maxSI_level; % evaluate linear function at maxSI_level

%% Plomp analysis
% A component
T_fit_emp_all_sii.A_pta = max(T_fit_emp_all_sii.pta_spl - srt_ref,0);

% D component
T_fit_emp_all_sii.D_srt_fsii = T_fit_emp_all_sii.fsii_indiv_srt - srt_ref - T_fit_emp_all_sii.A_pta;
T_fit_emp_all_sii.D_srt_fall = T_fit_emp_all_sii.fall_srt - srt_ref - T_fit_emp_all_sii.A_pta;
T_fit_emp_all_sii.D_srt_f1 = T_fit_emp_all_sii.f1_srt - srt_ref - T_fit_emp_all_sii.A_pta; % the same as SRT-PTA

%% Filter criteria
filter_crit_violet_all = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2 & T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c;

filter_crit_yellow_all = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI;

filter_crit_red_all = T_fit_emp_all_sii.case == 3 & T_fit_emp_all_sii.srt_pta_diff_80 >= thr_c;

%% Error analysis

% SRT based on 2 DP
p_srt2dp_lower = prctile(T_fit_emp_all_sii.srt_ci_2dp_lower(filter_crit_violet_all),0:25:100);
p_srt2dp_upper = prctile(T_fit_emp_all_sii.srt_ci_2dp_upper(filter_crit_violet_all),0:25:100);
p_srt2dp = (p_srt2dp_lower+p_srt2dp_upper)/2;
% average these two -> difference is only the SI error -> using the upper or lower distance to measured SI -> should be symmetric except for very high and low SI

% SRT based on 1 DP
p_srtsii_lower = prctile(T_fit_emp_all_sii.srt_ci_sii_lower(filter_crit_yellow_all),0:25:100);
p_srtsii_upper = prctile(T_fit_emp_all_sii.srt_ci_sii_upper(filter_crit_yellow_all),0:25:100);
p_srtsii = (p_srtsii_lower+p_srtsii_upper)/2;

% SRT based on NH fit
T_fit_emp_all_sii.f1_srt_range = T_fit_emp_all_sii.f1_srt - T_fit_emp_all_sii.f1_minSRT;

p_srtnh = prctile(T_fit_emp_all_sii.f1_srt_range(filter_crit_red_all),0:25:100);

% slope based on 2 DP
p_slp2dp_lower = prctile(T_fit_emp_all_sii.slope_error_2dp_lower(filter_crit_violet_all),0:25:100);
p_slp2dp_upper = prctile(T_fit_emp_all_sii.slope_error_2dp_upper(filter_crit_violet_all),0:25:100);
p_slp2dp = (p_slp2dp_lower+p_slp2dp_upper)/2;

p_slp2dp_rel = prctile(T_fit_emp_all_sii.slp_rel_error(filter_crit_violet_all),0:25:100);  % relative error, calculated on upper

% slope based on 1 DP
delta_s_h = sqrt(2)*sii_error;

% D component
T_fit_emp_all_sii.D_error_sii = T_fit_emp_all_sii.srt_ci_sii_upper - 5;
T_fit_emp_all_sii.D_error_2dp = T_fit_emp_all_sii.srt_ci_2dp_upper - 5;

p_Dsii = prctile(T_fit_emp_all_sii.D_error_sii(filter_crit_yellow_all),0:25:100);
p_D2dp = prctile(T_fit_emp_all_sii.D_error_2dp(filter_crit_violet_all),0:25:100);


%% Statistical analysis
% Comparison of distributions

% run statistical tests
% Mean: Welch test
% t-test, two-tail, unequal sample size, unequal variance
[h_mean_pta,p_mean_pta,ci_mean_pta,stats_mean_pta] = ttest2(T_fit_emp_all_sii.pta_spl(filter_crit_yellow_all),T_fit_emp_all_sii.pta_spl(filter_crit_violet_all),'alpha',0.05);

[h_mean_maxsi,p_mean_maxsi,ci_mean_maxsi,stats_mean_maxsi] = ttest2(T_fit_emp_all_sii.maxSI(filter_crit_yellow_all),T_fit_emp_all_sii.maxSI(filter_crit_violet_all),'alpha',0.05);

% Distribution -> Kolmogorov-Smirnoff test
% Kolmogorov-Smirnov-test, two-tail, unequal sample size
[h_dist_pta,p_dist_pta,stats_dist_pta] = kstest2(T_fit_emp_all_sii.pta_spl(filter_crit_yellow_all),T_fit_emp_all_sii.pta_spl(filter_crit_violet_all),'alpha',0.05);
[h_dist_maxsi,p_dist_maxsi,stats_dist_maxsi] = kstest2(T_fit_emp_all_sii.maxSI(filter_crit_yellow_all),T_fit_emp_all_sii.maxSI(filter_crit_violet_all),'alpha',0.05);

% Overlapping index -> see r2l_s6_plots_figA1.m

%% Prediction of SRT difference
rng(0)

T_fit_emp_all_sii.srt_diff = T_fit_emp_all_sii.fall_srt-T_fit_emp_all_sii.fsii_indiv_srt;

filter_crit = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2 & ~isnan(T_fit_emp_all_sii.srt_diff) & T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI;

T_pred_srt = T_fit_emp_all_sii(filter_crit,{'id_patient','sii_indiv_slope','sii_fit_dp_si2','srt_diff','fall_srt','fsii_indiv_srt'});
T_pred_srt.dsi_sign = T_pred_srt.sii_fit_dp_si2 - 50;

data_choice_prediction = {'sii_indiv_slope','dsi_sign'}; 
X = T_pred_srt{:,data_choice_prediction};
Y = T_pred_srt{:,{'srt_diff'}};


% partition data
K = 10;
cv = cvpartition(numel(Y), 'kfold',K);


mse = zeros(K,1);
for k=1:K
    % training/testing indices for this fold
    trainIdx = cv.training(k);
    testIdx = cv.test(k);

    % train GLM model
    mdl{k} = GeneralizedLinearModel.fit(X(trainIdx,:), Y(trainIdx), ...
        'linear');

    % prediction
    Y_hat = predict(mdl{k}, X(testIdx,:));
    Y_hat_all(testIdx) = Y_hat;

    % compute mean squared error
    mse(k) = mean((Y(testIdx) - Y_hat).^2);
end

% mean RMSE across folds
cv_rmse1 = mean(sqrt(mse));

% figure; plot(Y,Y_hat_all,'o') % predicted SRT diff

T_pred_srt.srt_diff_pred = Y_hat_all';
T_pred_srt.srt_sii_cor = T_pred_srt.fsii_indiv_srt + Y_hat_all';


T_fit_emp_all_sii = outerjoin(T_fit_emp_all_sii,T_pred_srt(:,{'id_patient','srt_diff_pred','srt_sii_cor'}),"Keys","id_patient","MergeKeys",true);

%% Evaluation
% arrange in table suitable for paper

T_glm = table();
for k = 1:K
    T_glm = [T_glm;array2table([ones(size(mdl{1}.Coefficients.pValue)),k*ones(size(mdl{1}.Coefficients.pValue)),[0:1:length(mdl{1}.Coefficients.pValue)-1]',mdl{k}.Coefficients.Estimate,mdl{k}.Coefficients.SE,mdl{k}.Coefficients.tStat,mdl{k}.Coefficients.pValue])];
    % k parameter estimate SE tStat pValue
end

% add averaged model parameters -> replace names in latex
mean_vec_est = [];
mean_vec_se = [];
mean_vec_tstat = [];
mean_vec_p = [];

for p = 1:length(mdl{1}.Coefficients.pValue)
    mean_vec_est = [mean_vec_est; mean(T_glm.Var4(T_glm.Var3==p-1 & T_glm.Var1==1))];
    mean_vec_se = [mean_vec_se; mean(T_glm.Var5(T_glm.Var3==p-1 & T_glm.Var1==1))];
    mean_vec_tstat = [mean_vec_tstat; mean(T_glm.Var6(T_glm.Var3==p-1 & T_glm.Var1==1))];
    mean_vec_p = [mean_vec_p; mean(T_glm.Var7(T_glm.Var3==p-1 & T_glm.Var1==1))];
end

T_glm = [T_glm; array2table([ones(size(mdl{1}.Coefficients.pValue)),(K+1)*ones(size(mdl{1}.Coefficients.pValue)),[0:1:length(mdl{1}.Coefficients.pValue)-1]',mean_vec_est,mean_vec_se,mean_vec_tstat,mean_vec_p])];

T_glm = renamevars(T_glm,["Var1","Var2","Var3","Var4","Var5","Var6","Var7"],["model","k","parameter","Estimate","SE","tStat","pValue"]);
T_glm.sig(T_glm.pValue<0.05) = 1;


% create latex table -> refine in tex doc
% sympref('FloatingPointOutput',1)
latex_table = latex(sym(round(T_glm{:,2:end},3)));


%% save data
writetable(T_fit_emp_all_sii,[res_folder filesep pres_paper filesep 's5_T_analysis.xlsx']);
save([res_folder filesep pres_paper filesep 's5_T_analysis.mat'],'T_fit_emp_all_sii');

save([res_folder filesep pres_paper filesep 'additional_variables_s5.mat'],'thr_c','srt_ref','sii_error','cv_rmse1');

