% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Local script 6: Plot figure 3  
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

filter_crit_red_all = (T_fit_emp_all_sii.case == 3 | T_fit_emp_all_sii.case == 5) & T_fit_emp_all_sii.srt_pta_diff_80 >= thr_c;  

%% Group data  

% C: Normal-hearing slope-based SRT estimation  
D_x = [...
       T_fit_emp_all_sii.pta_spl(filter_crit_red_all)];
D_y = [...
       T_fit_emp_all_sii.f1_srt(filter_crit_red_all)]- srt_ref;   
plot_colors = [...
               2*ones(size(T_fit_emp_all_sii.pta_spl(filter_crit_red_all)))];% ...
fill_colors = [...
               2*ones(size(T_fit_emp_all_sii.pta_spl(filter_crit_red_all)))];
 
idx_D = ones(size(D_x));

% markersize depending on frequency of same data point 
[D_unique,ia,ic] = unique([5*round(D_x/5),2*round(D_y/2)],'rows'); 
h_uni1 = hist(ic,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni1 = h_uni1/sum(h_uni1); 


% B: SII-based slope SRT estimation 
D_x = [T_fit_emp_all_sii.pta_spl(filter_crit_yellow_all)];
D_y = [T_fit_emp_all_sii.fsii_indiv_srt(filter_crit_yellow_all)]-srt_ref;
plot_colors = [3*ones(size(T_fit_emp_all_sii.pta_spl(filter_crit_yellow_all)))];
fill_colors = [3*ones(size(T_fit_emp_all_sii.pta_spl(filter_crit_yellow_all)))];
 
idx_D = ones(size(D_x));

% markersize depending on frequency of same data point 
[D_unique2,ia2,ic2] = unique([5*round(D_x/5),2*round(D_y/2)],'rows'); 
h_uni2 = hist(ic2,1:length(ia2)); % ia contains corresponding data points in the same order 
h_uni2 = h_uni2/sum(h_uni2); 


% A: Empirical slope-based SRT estimation 
D_x = [T_fit_emp_all_sii.pta_spl(filter_crit_violet_all)];
D_y = [T_fit_emp_all_sii.fall_srt(filter_crit_violet_all)]- srt_ref;
plot_colors = [4*ones(size(T_fit_emp_all_sii.pta_spl(filter_crit_violet_all)))];
fill_colors = [4*ones(size(T_fit_emp_all_sii.pta_spl(filter_crit_violet_all)))];
 
idx_D = ones(size(D_x));

% markersize depending on frequency of same data point 
[D_unique3,ia3,ic3] = unique([5*round(D_x/5),2*round(D_y/2)],'rows'); 
h_uni3 = hist(ic3,1:length(ia3)); % ia contains corresponding data points in the same order 
h_uni3 = h_uni3/sum(h_uni3); 

%% Figure  

figh2 = figure; 

subplot(1,3,1)
hold on;
box on;

line([-5 120],[-5 120]-srt_ref,[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

scatter(D_unique3(:,1),D_unique3(:,2),4000*log10(h_uni3+1)','MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:)); 

xlabel('PTA [dB SPL]')
ylabel('SRT_f - SRT_{NH} [dB]') 
title('Empirical slope')
axis([-5 120 0 80])
text(0,76,['{\itN} = ' num2str(length(ic3))])

%
subplot(1,3,2)
hold on;
box on;

line([-5 120],[-5 120]-srt_ref,[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

scatter(D_unique2(:,1),D_unique2(:,2),4000*log10(h_uni2+1)','MarkerEdgeColor',fig_properties.fillcolors(3,:),'MarkerFaceColor',fig_properties.linecolors2(3,:)); 
xlabel('PTA [dB SPL]')
ylabel('SRT_h - SRT_{NH} [dB]') 
title('SII-based slope')
axis([-5 120 0 80])
text(0,76,['{\itN} = ' num2str(length(ic2))])
%
subplot(1,3,3)

hold on;
box on;

line([-5 120],[-5 120]-srt_ref,[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

scatter(D_unique(:,1),D_unique(:,2),4000*log10(h_uni1+1)','MarkerEdgeColor',fig_properties.fillcolors(2,:),'MarkerFaceColor',fig_properties.linecolors2(2,:)); 

xlabel('PTA [dB SPL]')
ylabel('SRT_n - SRT_{NH} [dB]')  
title('Normal-hearing slope')
axis([-5 120 0 80])
text(0,76,['{\itN} = ' num2str(length(ic))])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.5*figw figh]);

if sflag
    print(figh2,[fig_folder filesep pfig_paper filesep 'Fig3_srt_pta_subplot.eps'],'-painters','-depsc','-r300');
    print(figh2,[fig_folder filesep pfig_paper filesep 'Fig3_srt_pta_subplot.png'],'-dpng','-r300');
    
end

