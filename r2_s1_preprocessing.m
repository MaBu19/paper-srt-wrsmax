% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Remote script 1: preprocessing 
%
% Mareike Buhl
% mareike.buhl@pasteur.fr
%

clear all;
close all;
clc;

% change here to have all paths correct:
dataflag = 1; % 1: remote analysis (EK), 2: synthetic example analysis (MB)

% if complete data is not available but results should be loaded 
readflag = 1; % 1: default, 0: disable data extraction

warning('off');

%% read data
if dataflag == 1
    if readflag == 1
        % paths
        addpath(genpath('O:\Github\h-ol-remote-example-analysis\functions_remote\'));
        data_folder = 'O:\Projekte\bigData\BigData_Oldenburg\data\';
        res_folder = 'O:\Github\h-ol-remote-example-analysis\results_remote\';
        fig_folder = 'O:\Github\h-ol-remote-example-analysis\figures_remote';

        % read data
        % Make connection to database
        conn = database('CASwarehouse','','');

        % Set query to execute on the database
        query = ['SELECT * FROM CASwarehouse.dbo.patient_audiogram'];

        % Execute query and fetch results
        t_pat_ag = fetch(conn,query);

        % Close connection to database
        close(conn)

        % Clear variables
        clear conn query


    else
        addpath(genpath('./functions_remote/'));
        res_folder = './results_remote/';
        fig_folder = './figures_remote/';
    end
elseif dataflag == 2
    % paths
    path_str = '_remote'; % '_remote'|''
    addpath(genpath(['./functions' path_str '/']));
    data_folder = '../../projects/hanover-remote-analysis/h-ol-example-data/CASwarehouse_example/';

    res_folder = './results_example/';
    fig_folder = './figures_example/';

    % read data
    f1 = dir(fullfile(data_folder,'patient_audiogram.mat')); 
    t_patient_audiogram = load([data_folder f1.name]);
    t_pat_ag = t_patient_audiogram.patient_audiogram;

end

pfig = 'paper'; % subfolder for this script

if ~exist([res_folder filesep pfig], 'dir')
    mkdir([res_folder filesep pfig]);
end

if ~exist([fig_folder filesep pfig], 'dir')
    mkdir([fig_folder filesep pfig]);
end

T_data = renamevars(t_pat_ag,["AC_monosyllabic_60dB","AC_monosyllabic_80dB","AC_monosyllabic_100dB","AC_monosyllabic_110dB"],["x60","x80","x100","x110"]);

%% Delete patients' ears that have no audiogram or no speech test 

N_all = length(unique(T_data.id_patient)); % number of patients before any preprocessing step 

% number of data points FMST  
si_vec = T_data{:,{'x60','x80','x100','x110'}};  
idxy = ~isnan(si_vec);
T_data.m_dp = sum(idxy,2);

% number of data points audiogram 
ag_vec = T_data{:,["AC_250_Hz","AC_500_Hz","AC_1000_Hz","AC_1500_Hz","AC_2000_Hz","AC_3000_Hz","AC_4000_Hz","AC_6000_Hz","AC_8000_Hz"]};
idxa = ~isnan(ag_vec);
T_data.m_ag = sum(idxa,2);

% save histograms of data points before removing the missing data -> to
% mention numbers in paper if interesting 
m_dp_histo_pre = hist(T_data.m_dp,[0 1 2 3]);
m_ag_histo_pre = hist(T_data.m_ag,0:9);

T_data(T_data.m_dp==0,:) = []; 
T_data(T_data.m_ag==0,:) = []; 

 
%% Interpolate audiograms 
 
T_data = ag_interpolate(T_data,{'AC_250_Hz','AC_500_Hz','AC_1000_Hz','AC_1500_Hz','AC_2000_Hz','AC_3000_Hz','AC_4000_Hz','AC_6000_Hz','AC_8000_Hz'}); 

%% Calculate PTA    
T_data.pta = nanmean(T_data{:,["AC_500_Hz","AC_1000_Hz","AC_2000_Hz","AC_4000_Hz"]},2);
T_data.pta_bc = nanmean(T_data{:,["BC_500_Hz","BC_1000_Hz","BC_2000_Hz","BC_4000_Hz"]},2);
T_data.pta_abg = T_data.pta-T_data.pta_bc; 
T_data.pta_spl = nanmean(db_hl_spl([500 1000 2000 4000],T_data{:,["AC_500_Hz","AC_1000_Hz","AC_2000_Hz","AC_4000_Hz"]},'HL'),2); 

T_data.A_05 = db_hl_spl(500,T_data{:,["AC_500_Hz"]},'HL'); 
T_data.A_1 = db_hl_spl(1000,T_data{:,["AC_1000_Hz"]},'HL');
T_data.A_2 = db_hl_spl(2000,T_data{:,["AC_2000_Hz"]},'HL');
T_data.A_4 = db_hl_spl(4000,T_data{:,["AC_4000_Hz"]},'HL'); 

  
%% Keep the better ear of each patient
[g_pat,id] = findgroups(T_data.id_patient);
pta_better = splitapply(@min,T_data.pta,g_pat); % max just takes one of two equal values
t_ag_filter = table(id,pta_better,'VariableNames',{'id_patient','pta'});
T_data = innerjoin(T_data,t_ag_filter,"Keys",{'id_patient','pta'});  
% also cut some last double ears (not removed due to same pta)
[~,pa,~] = unique(T_data.id_patient,'first');   
T_data = T_data(pa,:);


%% Calculate Bisgaard class

ag_freq = [250 500 1000 1500 2000 3000 4000 6000 8000]; 
for n = 1:height(T_data)
    [bisg_calc(n), ag_rmse(n)] = calc_bisg_rmse(T_data{n,["AC_250_Hz","AC_500_Hz","AC_1000_Hz","AC_1500_Hz","AC_2000_Hz","AC_3000_Hz","AC_4000_Hz","AC_6000_Hz","AC_8000_Hz"]}, ag_freq);
end
T_data.bisg = bisg_calc'; 
T_data.bisg_rmse = ag_rmse'; 
T_data.freqdiff = nanmean([T_data.AC_2000_Hz,T_data.AC_4000_Hz],2)-nanmean([T_data.AC_500_Hz,T_data.AC_1000_Hz],2);  
  
% figure; 
% plot(T_data{T_data.bisg == "S1" & T_data.bisg_rmse < 15,["AC_250_Hz","AC_500_Hz","AC_1000_Hz","AC_1500_Hz","AC_2000_Hz","AC_3000_Hz","AC_4000_Hz","AC_6000_Hz","AC_8000_Hz"]}'); axis ij 
% % TO DO: inspect later: T_data.bisg_rmse > 25 (can still be excluded later
% % in analysis based on error) 
% 
% figure; 
% hist(T_data.bisg,unique(bisg_calc))
% 
% figure; 
% hist(T_data.bisg_rmse,0:5:50)

%% additional histograms -> after filtering the data used for analysis 
channel_histo = hist(categorical(T_data.channel),{'L','R'}); 

gender_histo = hist(categorical(T_data.gender),{'F','M'}); 

age_histo = hist(T_data.age_measurement,0:5:120); 
age_percentiles = prctile(T_data.age_measurement,0:5:100); 

T_data.test_date = datetime(T_data.test_date);
T_data.year = year(T_data.test_date); 
date_histo = hist(T_data.year,min(T_data.year):1:max(T_data.year)); 

min_year = min(T_data.year); 
max_year = max(T_data.year); 

%% Output table for SII calculation 
% -> only unique audiograms to save computing time + correspondance table 

[~,id_unique,idx_back] = unique(T_data{:,["AC_250_Hz","AC_500_Hz","AC_1000_Hz","AC_1500_Hz","AC_2000_Hz","AC_3000_Hz","AC_4000_Hz","AC_6000_Hz","AC_8000_Hz"]},'rows'); 
% later use to project back SII results to patients with same audiogram

T_AG_SII = T_data(id_unique,["id_patient","AC_250_Hz","AC_500_Hz","AC_1000_Hz","AC_1500_Hz","AC_2000_Hz","AC_3000_Hz","AC_4000_Hz","AC_6000_Hz","AC_8000_Hz"]);

T_data.idx_back = idx_back; 
% T_AG_SII.idx_back = idx_back(id_unique); % would also be correct
% correspondance, but neither needed nor used in SII script 
% -> in r2_s4_srt_assumptions_v3: join tables (data and SII slopes) via idx_back

%% Output table for SRT estimation 
T_all = T_data(:,{'id_patient','pta','pta_spl','pta_abg','m_ag','A_05','A_1','A_2','A_4','bisg','bisg_rmse','freqdiff','m_dp','x60','x80','x100','x110','idx_back'}); 

%% Save 
% for SII calculations 
writetable(T_AG_SII,[res_folder filesep pfig filesep 's1_T_AG_SII.xlsx']);
save([res_folder filesep pfig filesep 's1_T_AG_SII.mat'],'T_AG_SII','id_unique','idx_back');

% for complete analysis (script s4)
writetable(T_all,[res_folder filesep pfig filesep 's1_T_all.xlsx']);
save([res_folder filesep pfig filesep 's1_T_all.mat'],'T_all');

% additional variables preprocessing 
save([res_folder filesep pfig filesep 's1_additional_variables.mat'],'N_all','min_year','max_year','m_dp_histo_pre','m_ag_histo_pre','channel_histo','gender_histo','age_histo','age_percentiles','date_histo');


