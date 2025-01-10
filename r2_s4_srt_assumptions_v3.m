% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Remote script 4: SRT estimation  
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

pfig = 'paper'; % subfolder for this script
if ~exist([res_folder filesep pfig], 'dir')
    mkdir([res_folder filesep pfig]);
end
if ~exist([fig_folder filesep pfig], 'dir')
    mkdir([fig_folder filesep pfig]);
end


%% load T_all from s1
load([res_folder pfig filesep 's1_T_all']);

%% Parameters
SRTref_FMSTq = 29.3; % ref SRT quiet

% load SII slope for individual audiogram 
T_model_indiv = readtable([res_folder pfig filesep 's3_T_slope_indiv_sii-hl0-LD-mhh_opt_v2.txt']);

T_model_indiv = renamevars(T_model_indiv,{'id_patient','slope'},{'bisgaard','SII_slope'});
T_model_indiv.rel_slope = T_model_indiv.SII_slope/0.0307; % NH SII slope: 0.0307
for n = 1:height(T_model_indiv)
    a = char(T_model_indiv.bisgaard(n));
    T_model_indiv.id_patient(n) = str2num(a(2:end));
end

% add slope to T_all 
T_model_indiv = join(T_model_indiv,T_all(:,{'id_patient','idx_back'}),'Keys','id_patient');
T_all = innerjoin(T_all,T_model_indiv(:,{'idx_back','SII_slope','r2','rel_slope'}),'Keys','idx_back');
T_all.sii_indiv_slope = T_all.rel_slope*4.5;

% SI confidence intervals Holube et al. (2018)
% (95% confidence interval -> significant difference expectable if two
% results differ by the given amount) , table 5, one test list
% - si_vector 0:5:100 -> index: SI/5+1
lower_ci = [0 1.1 3.4 6.2 9.4 12.8 16.5 20.3 24.4 28.6 32.9 37.4 42 46.8 51.8 57 62.4 68 74 80.6 88.3];
upper_ci = [11.7 19.4 26 32 37.6 43 48.2 53.2 58 62.6 67.1 71.4 75.6 79.7 83.5 87.2 90.6 93.8 96.6 98.9 100];

% SII error 
sii_error = 0.00084;  % estimated with 10 SII runs for Bisgaard audiograms and different levels (max of std of all conditions) 


%% SI calculations
% calculate maximum discrimination for every patient
level_vec = [60 80 100 110];
[SI_max,level_max,SI_min,level_min] = max_SI_discr(T_all{:,{'x60','x80','x100','x110'}},level_vec);
[SI_level_max, SI_level_max_ia, SI_level_max_ic]  = unique([SI_max level_max'],'rows');

T_all.maxSI = SI_max;
T_all.maxSI_level = level_max';
T_all.maxSI_level_combi = SI_level_max_ic;

T_all.minSI = SI_min;
T_all.minSI_level = level_min';

%% Fit linear function to SI results in constant slope area

% level_vec: taken from above - property of data set
level_combis = nan(height(T_all),length(level_vec)); % to store the available level combinations for all psychometric functions
si_combis = nan(height(T_all),length(level_vec));
for n = 1:height(T_all)
    
    if mod(n,1000) == 0
        disp(['Progress: ' num2str(round(n/height(T_all)*100,1)) ' %']);
    end

    % sort input for fit
    si_vec = T_all{n,{'x60','x80','x100','x110'}}';
    idxy = ~isnan(si_vec);
    level_combis(n,:) = idxy;
    si_combis(n,:) = si_vec';

    % special cases of data availability -> how to fix fit for these cases:
    % 1) cut si_vec after maximum (to exclude decreasing part if present)
    [valmax,idxmax]=max(si_vec);
    si_vec(idxmax:end) = valmax; % replace decreasing part by maximum - only for fit, to guide fit to maxSI
    si_vec(~idxy) = nan; % points which were nan before should stay nan

    % additional restriction to datapoints at constant part of assumed
    % psychometric function
    SI_limit = 0.85*T_all.maxSI(n);
    idx_limit = si_vec<=SI_limit;
    % + same with lower limit
    SI_limit2 = 0.15*T_all.maxSI(n);
    idx_limit2 = si_vec>=SI_limit2;

    si_vec(~idx_limit) = nan;
    si_vec(~idx_limit2) = nan;

    T_all.m_dp_fit(n) = sum(idx_limit & idx_limit2);

    levels_fit = level_vec(idx_limit & idx_limit2)';
    si_fit = si_vec(idx_limit & idx_limit2);

    % choose one of the two datapoints (audible and closer to 50% SI)
    idx_aud = levels_fit > T_all.pta_spl(n)-10; % 10 dB difference: assumption that one can hear a bit below the PTA, also due to generally a bit sloping hearing losses

    [mlev,idx_m2] = min(abs(si_fit(idx_aud)-50));  % index used for fit  
    if length(si_fit)>1
        if idx_aud(1:2) == [0 1]'
            idx_m2 = idx_m2+1; % in this case idx_m2 = 1 before, but idx should be 2 since only second data point is audible
        end

    end

    if ~isempty(idx_m2)
        T_all.sii_fit_dp_lev2(n) = levels_fit(idx_m2);
        T_all.sii_fit_dp_si2(n) = si_fit(idx_m2);
        T_all.sii_fit_dp_mlev(n) = mlev;
    else
        T_all.sii_fit_dp_lev2(n) = nan;
        T_all.sii_fit_dp_si2(n) = nan;
        T_all.sii_fit_dp_mlev(n) = nan;
    end


    % calculate errors for different patient cases -> related to SI according
    % to Holube; difference between measured SI and values given above in lower_ci and upper_ci
    % case 2 -> one dp as defined by idx_m
    % case 2, 2 DP -> two datapoints as given in si_fit -> sufficient to save these two errors and idx_m since case 3 datapoint is a subset of the two dp here
    si_errors_lower = abs(lower_ci(si_fit/5+1)-si_fit'); % upper und lower refers to error bounderies, for both data points, respectively
    si_errors_upper = abs(upper_ci(si_fit/5+1)-si_fit');

    %% Case x1: 

    % linear fit to all data points (2) in the chosen SI range -> slope is variable (empirical slope based SRT estimation)
    [ict_fitted(n,1),slp_fitted(n,1),srt_fitted(n,1),gof_tmp] = fit_linear_slope(level_vec(idx_limit & idx_limit2)',si_vec(idx_limit & idx_limit2),1);
    T_all.fall_gof(n) = gof_tmp;
    T_all.fall_sse(n) = gof_tmp.sse;

    % slope error for fit with two dp
    if length(si_fit)>1
        T_all.slope_error_2dp_lower(n) = 1/(levels_fit(2)-levels_fit(1))*(si_errors_lower(2)+si_errors_lower(1));
        T_all.slope_error_2dp_upper(n) = 1/(levels_fit(2)-levels_fit(1))*(si_errors_upper(2)+si_errors_upper(1));

    else
        T_all.slope_error_2dp_lower(n) = nan;
        T_all.slope_error_2dp_upper(n) = nan;
    end

    % interesting to calculate relative error in addition: 
    T_all.slp_rel_error(n,1) = T_all.slope_error_2dp_upper(n)/slp_fitted(n,1); 

    % srt error based on slope error estimated above -> note:
    % same data point chosen for SI error as for SII-fit 
    if length(si_fit)>1 && ~isempty(idx_m2)
        T_all.srt_ci_2dp_lower(n) = 1/slp_fitted(n,1)*si_errors_lower(idx_m2) + (abs(si_fit(idx_m2)-50))*(1/slp_fitted(n,1)^2)*T_all.slope_error_2dp_lower(n); % before idx_m
        T_all.srt_ci_2dp_upper(n) = 1/slp_fitted(n,1)*si_errors_upper(idx_m2) + (abs(si_fit(idx_m2)-50))*(1/slp_fitted(n,1)^2)*T_all.slope_error_2dp_upper(n);
    else
        T_all.srt_ci_2dp_lower(n) = nan;
        T_all.srt_ci_2dp_upper(n) = nan;
    end

    %% Case x2:
    % linear fit with assumed slope as relative to SII slope for individual
    % audiogram (SII-based slope SRT estimation)

    % srt error for fit with one dp and SII-indiv-slope -> same data point as below for fit
    if ~isempty(idx_m2)
        T_all.srt_ci_sii_lower(n) = 1/T_all.sii_indiv_slope(n)*si_errors_lower(idx_m2) + (abs(si_fit(idx_m2)-50))*(1/T_all.sii_indiv_slope(n)^2)*sii_error;
        T_all.srt_ci_sii_upper(n) = 1/T_all.sii_indiv_slope(n)*si_errors_upper(idx_m2) + (abs(si_fit(idx_m2)-50))*(1/T_all.sii_indiv_slope(n)^2)*sii_error;
        % (dominated by SI error)
    else
        T_all.srt_ci_sii_lower(n) = nan;
        T_all.srt_ci_sii_upper(n) = nan;
    end

    % fit
    [ict_fitted(n,2),slp_fitted(n,2),srt_fitted(n,2),gof_tmp] = fit_linear_slope(levels_fit(idx_m2),si_fit(idx_m2),2,T_all.sii_indiv_slope(n));

    T_all.fsii_indiv_gof(n) = gof_tmp;
    T_all.fsii_indiv_sse(n) = gof_tmp.sse;

    %% Case x3: 
    % simplest assumption to fit psychometric function with normal-hearing slope and maxSI = 100 % - if nothing else is possible given the available data

    [srt_psyfun(n,1),~,maxsi_psyfun(n,1),gof_tmp] = fit_srt_psyfun(T_all.maxSI_level(n),T_all.maxSI(n),1);
    T_all.f1_gof(n) = gof_tmp;
    T_all.f1_sse(n) = gof_tmp.sse;

    % Error: Best we can do: estimate range of potential SRT -> since real slope is assumed to be lower than NH slope, the real SRT can only be lower than the estimated SRT with the NH fit; lower bound is NH SRT or audibility (here same assumption as used later: SRTs are kept until 10 dB lower than pta, otherwise rejected due to audibility problem)
    T_all.f1_minSRT(n) = max(SRTref_FMSTq,T_all.pta_spl(n)-10); % range calculated in local script

end

% assemble results to table
T_all.fall_intercept = ict_fitted(:,1);
T_all.fall_slope = slp_fitted(:,1);
T_all.fall_srt = srt_fitted(:,1);

T_all.fsii_indiv_intercept = ict_fitted(:,2);
T_all.fsii_indiv_slope = slp_fitted(:,2);
T_all.fsii_indiv_srt = srt_fitted(:,2);

% 1 parameter psychometric function fit (case 3)
T_all.f1_srt = srt_psyfun(:,1);
T_all.f1_slope = 4.5*ones(size(T_all.f1_srt));
T_all.f1_maxSI = maxsi_psyfun(:,1);



%% add additional variables to table
% available level combinations across patients
sum(level_combis);
[level_combis_all,~,ic] = unique(level_combis,'rows');
% add index for level combination to fit result table
T_all.level_combi = ic;

si_combis(isnan(si_combis)) = 222; % to be able to apply unique
[si_combis_all,~,is] = unique(si_combis,'rows');
T_all.si_combi = is;

% indices for data availability/patient groups for which different SRT
% estimation or assumptions are possible (cases)
T_all.c_fall_srt(~isnan(T_all.fall_srt)) = 1;
T_all.c_fsii_indiv_srt(~isnan(T_all.fsii_indiv_srt)) = 1;

% case numbers:
T_all.case(T_all.c_fall_srt == 1) = 1; % SRT can be estimated - two datapoints available for fit (subset of SII-fit)
T_all.case(T_all.c_fsii_indiv_srt == 1) = 2; % SRT can be estimated - one datapoint available for fit - only fit with SII-based slope possible
T_all.case(T_all.m_dp == 1 & (T_all.maxSI >= 80)) = 3; % SRT cannot be estimated since only one datapoint is available, corresponding to maxSI already
T_all.case(T_all.maxSI < 50*10/9) = 4; % if maxSI too small, not important how many datapoints were measured -> no SRT since maxSI too low
T_all.case(T_all.m_dp == 1 & (T_all.maxSI < 80) & (T_all.maxSI >= 50*10/9)) = 5; % 1 DP and lower maxSI (between 55 and 80%) -> also patients for which 1 DP fit can be executed (together with case 3, with the difference that case 3 fulfils hearing aid indication criteria)
T_all.case(T_all.c_fsii_indiv_srt == 0 & T_all.m_dp_fit == 1 & (T_all.maxSI >= 50*10/9)) = 6; % for these patients one datapoint was available for fit, but due to the low SII slope (and in combination with the measured level and obtained SI), 50% was not achieved and therefore no SRT was obtained -> only possible way to estimate SRT is based on NH slope (which potentially overestimates the real SRT if the real slope is lower)

%% individual psychometric function as checkplot - some random examples

T_plot = T_all(T_all.case == 2 & T_all.m_dp_fit == 2,:);
% T_plot = T_all(T_all.case == 3 ,:);
xvec = 10:110;

if height(T_plot) > 0
for n = 1 %:height(T_plot)
    if T_plot.m_dp_fit(n) >= 1

        figh1 = figure;
        title([num2str(T_plot.id_patient(n)) ', ' T_plot.bisg(n)] )
        hold on;
        box on;
        % normative data
        plot(xvec,psyfun(29.3,4.5/100,xvec,100),'linewidth',1.5)

        plot(xvec,T_all.fall_intercept(T_all.id_patient == T_plot.id_patient(n))+T_all.fall_slope(T_all.id_patient == T_plot.id_patient(n))*xvec,'linewidth',1.5)
        plot(xvec,T_all.fsii_indiv_intercept(T_all.id_patient == T_plot.id_patient(n))+T_all.fsii_indiv_slope(T_all.id_patient == T_plot.id_patient(n))*xvec,'linewidth',1.5)

        % add PTA to plot
        line([T_plot.pta_spl(n) T_plot.pta_spl(n)],[0 100],[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

        % SRT
        plot(T_plot.fall_srt(n),50,'s','linewidth',1.5)
        plot(T_plot.fsii_indiv_srt(n),50,'d','linewidth',1.5)

        % empirical data
        plot(level_vec,T_all{T_all.id_patient==T_plot.id_patient(n),{'x60','x80','x100','x110'}},'o','linewidth',1.5)

        % plot properties
        xlabel('Speech level [dB SPL]')
        ylabel('Speech intelligibility [%]')
        axis([10 110 0 100])

        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 12*0.85 11*0.85]);

        if 1
            print(figh1,[fig_folder filesep pfig filesep 'remote_psyfun-ex-' num2str(T_plot.id_patient(n)) '.png'],'-dpng','-r300');
        end
    end
end
end

%% save tables (without individual raw data)

T_fit = removevars(T_all,{'x60','x80','x100','x110'}); % remove empirical data 
writetable(T_fit,[res_folder filesep pfig filesep 's4_T_fit.xlsx']); 
save([res_folder filesep pfig filesep 's4_T_fit.mat'],'T_fit');
