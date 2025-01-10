% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Local script 6: Plot figure 7 and 8 (appendix) 
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
pfig_paper = 'paper_plots_main';

if ~exist([fig_folder filesep pfig_paper], 'dir')
    mkdir([fig_folder filesep pfig_paper]);
end

sflag = 1;
 
%% load data (remote analysis)

load([res_folder pres_paper filesep 's5_T_analysis']);  
load([res_folder pres_paper filesep 'additional_variables_s5']) 

%% figure properties 
figw = 12*0.85; 
figh = 11*0.85; 

map2 = [colormap('lines'); ... 
        1 1 1; ... % white - 257
        1 0.55 0]; % orange - 258

fig_properties.linecolors2 = map2;  
fig_properties.fillcolors = map2;
fig_properties.fillcolors(2,:) = [255 114 118]/255; 
fig_properties.fillcolors(3,:) = [241 235 156]/255; 
fig_properties.fillcolors(4,:) = [177 162 250]/255;  


%% filter criteria 
filter_crit_violet = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2 & ~isnan(T_fit_emp_all_sii.srt_diff) & T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI; 
filter_crit_violet_all = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2 & T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c; 

filter_crit_yellow_all = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI; 

filter_crit_red_all = T_fit_emp_all_sii.case == 3 & T_fit_emp_all_sii.srt_pta_diff_80 >= thr_c; 

%% Figure 7: PTA distributions 

T_stats = table(); 

pta_bins = -10:5:120; 

figh5 = figure('Visible','on'); 

hold on; box on; 
h3 = hist(T_fit_emp_all_sii.pta_spl(filter_crit_yellow_all),pta_bins);  
h4 = hist(T_fit_emp_all_sii.pta_spl(filter_crit_violet_all),pta_bins); 

subplot(2,1,1)
b1 = bar(pta_bins,h3'); 
xlim([pta_bins(1)-3 pta_bins(end)+3])
legend({'Half-determined'},'Location','NorthWest')  
title('PTA for fully-determined and half-determined patients')
ylabel('Number of patients')

subplot(2,1,2)
b2 = bar(pta_bins,h4'); 
xlim([pta_bins(1)-3 pta_bins(end)+3])
legend({'Fully-determined'},'Location','NorthWest')  

b1.FaceColor = map2(3,:);
b2.FaceColor = map2(4,:);

xlabel('PTA [dB SPL]')
ylabel('Number of patients')
xlim([pta_bins(1)-3 pta_bins(end)+3])

text(-0.3,max(sum([h3; h4]))*0.97,['{\itN} = ' num2str(height(T_fit_emp_all_sii(filter_crit_yellow_all,:)))])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figw figh]);

if sflag
    print(figh5,[fig_folder filesep pfig_paper filesep 'FigA1A_histo_PTA_yellow_violet.eps'],'-painters','-depsc','-r300');
    print(figh5,[fig_folder filesep pfig_paper filesep 'FigA1A_histo_PTA_yellow_violet.png'],'-dpng','-r300');
    
end

% Overlapping index 
T_stats.variable(1) = categorical({'PTA'}); 
T_stats.eta_yv(1) = sum(min(h3/sum(h3),h4/sum(h4))); 

%% Figure 8: WRS_max distributions 

maxsi_bins = 0:5:100; 

figh5 = figure('Visible','on'); 

hold on; box on; 
h3 = hist(T_fit_emp_all_sii.maxSI(filter_crit_yellow_all),maxsi_bins);
h4 = hist(T_fit_emp_all_sii.maxSI(filter_crit_violet_all),maxsi_bins);  

subplot(2,1,1)
b1 = bar(maxsi_bins,h3'); 
xlim([50 103])

legend({'Half-determined'},'Location','NorthWest') 

title('WRS_{max} for fully-determined and half-determined patients')
ylabel('Number of patients')

subplot(2,1,2)
b2 = bar(maxsi_bins,h4'); 
xlim([maxsi_bins(1)-3 maxsi_bins(end)+3])
legend({'Fully-determined'},'Location','NorthWest')  


b1.FaceColor = map2(3,:);
b2.FaceColor = map2(4,:);

xlabel('WRS_{max} [%]')
ylabel('Number of patients')
xlim([50 103])
text(-0.3,max(sum([h3; h4]))*0.97,['{\itN} = ' num2str(height(T_fit_emp_all_sii(T_fit_emp_all_sii.case == 2,:)))])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figw figh]);

if sflag
    print(figh5,[fig_folder filesep pfig_paper filesep 'FigA1B_histo_maxSI_yellow_violet.eps'],'-painters','-depsc','-r300');
    print(figh5,[fig_folder filesep pfig_paper filesep 'FigA1B_histo_maxSI_yellow_violet.png'],'-dpng','-r300');
    
end

% Overlapping index 
T_stats.variable(2) = 'maxSI'; 
T_stats.eta_yv(2) = sum(min(h3/sum(h3),h4/sum(h4))); 





