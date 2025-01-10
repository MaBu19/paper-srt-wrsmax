% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Remote script 3: Combination of SII results to one table  
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

pfig = 'paper'; 
res_combi_folder = [res_folder filesep pfig]; % subfolder for this script
if ~exist(res_combi_folder, 'dir')
    mkdir(res_combi_folder);
end

res_folders = { [res_folder pfig '/s2_sii/hl0-LD-mhh_opt_v2']}; 
agtypes = {'_indiv'}; 

% define empty table: all columns as in T_olsa + folder name
for agtype = 1:length(agtypes)
        fname_part = [agtypes{agtype}]; 

T_all = table;

% list of folders which include data to be included in the table 
for k = 1:length(res_folders)
 
    fname_part_tmp = strrep(strrep(res_folders{k},[res_folder pfig '/s2_'],''),'/','-'); 
    fname_part = [fname_part, '_', fname_part_tmp];

    % per folder, for all files matching the file name criteria:
    % - read data
    % - split filename -> info should be stored in table as well
    % - add all data to table

    files = dir(fullfile([res_folders{k} ],['sii_*' agtypes{agtype} '.txt'])); 

    for f = 1:length(files)
        T_data = readtable([res_folders{k} '/' files(f).name ]);
        st_prop = strsplit(files(f).name,'.');
        st_str = strsplit(st_prop{1},'_');
        T_prop = array2table(st_str);
        T_folder = table; T_folder.folder = res_folders(k);
        T_line = [T_folder,T_prop,T_data];
        T_all = [T_all;T_line];

    end 
end

T_all = renamevars(T_all,["st_str1","st_str2","st_str3","st_str4","st_str5","st_str6","st_str7","st_str8","st_str9"],["model","speechtest","bisgaard","language","speechlevel","noisetype","noiselevel","q_n","type"]);


% filename_out: includes the list of folders
writetable(T_all,[res_combi_folder '/s3_T_summary' fname_part],FileType='text')
end 

%% do the same for slope tables: 
% define empty table: all columns as in T_slope + folder name
for agtype = 1:length(agtypes)
        fname_part = [agtypes{agtype}]; 

T_all_slope = table;

% list of folders which include data to be included in the table -> loop
for k = 1:length(res_folders)

    fname_part_tmp = strrep(strrep(res_folders{k},[res_folder pfig '/s2_'],''),'/','-'); 
    fname_part = [fname_part, '_', fname_part_tmp];

    % per folder, for all files matching the file name criteria:
    % - read data
    % - split filename -> info should be stored in table as well
    % - add all data to table

    files = dir(fullfile([res_folders{k} ],['Tslope_*' agtypes{agtype} '_p*.txt'])); %Tslope_sii_FMST_DE_ccitt_-50_quiet_indiv_p3505

    for f = 1:length(files)
        T_data = readtable([res_folders{k} '/' files(f).name ]);
        T_all_slope = [T_all_slope;T_data];

    end
end

% filename_out: includes the list of folders
writetable(T_all_slope,[res_combi_folder '/s3_T_slope' fname_part],FileType='text')
end 


%% check: SII error (dmv) varying across Bisgaard class or speech level?
% okay to assume one value for all conditions? 

% mean(T_all.dmv)
% std(T_all.dmv)
% max(T_all.dmv) % 0.00084 -> used for all conditions in r2_s4_srt_assumptions_v3 

