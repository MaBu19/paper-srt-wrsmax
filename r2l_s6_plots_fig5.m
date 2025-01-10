% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Local script 6: Plot figure 5  
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
fig_properties.fillcolors(4,:) = [177 162 250]/255; 

%% filter criteria 
filter_crit_violet = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2 & ~isnan(T_fit_emp_all_sii.srt_diff) & T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI; 

%% Figure 5: Prediction of SRT difference between the two SRT estimation procedures 
% first model
[P1,xbins] = calc_percentiles_sections(T_fit_emp_all_sii.srt_diff,T_fit_emp_all_sii.srt_diff_pred,-15:2.5:15,filter_crit_violet,10:20:90); 

D_x = T_fit_emp_all_sii.srt_diff(filter_crit_violet); 
D_y = T_fit_emp_all_sii.srt_diff_pred(filter_crit_violet); 

plot_colors = 4*ones(size(D_x)); 

idx_D = 1:length(D_x); 


% markersize depending on frequency of same data point 
[D_unique,ia,ic] = unique([1*round(D_x/1),1*round(D_y/1)],'rows'); 
h_uni1 = hist(ic,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni1 = h_uni1/sum(h_uni1); 

[r2,p,bias,rmsv,rl2,ru2] = statistical_analysis(D_x(idx_D),D_y(idx_D));

figh6 = figure('Visible','on');
hold on;
box on;

line([-15 120],[-15 120],[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

% "legend"
legend_size1 = [round(max(h_uni1)*length(ic)) round(max(h_uni1)*length(ic)/2) round(max(h_uni1)*length(ic)/4) round(max(h_uni1)*length(ic)/8) round(max(h_uni1)*length(ic)/16)]; 
legend_y1 = -5; 
for n = 1:length(legend_size1)
    scatter(10,legend_y1-(n-1)*3,4000*log10(legend_size1(n)/length(ic)+1),'MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:));  
    text(12,legend_y1-(n-1)*3,num2str(round(legend_size1(n))))
end 

lwidth = [0.5 1.5 2.5 1.5 0.5];
for n = 1:size(P1,2)
    plot(xbins,P1(:,n),'-','color',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.fillcolors(4,:),'linewidth',lwidth(n))
end

scatter(D_unique(:,1),D_unique(:,2),4000*log10(h_uni1+1)','MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:)); 

xlabel('SRT difference (observed) [dB]')
ylabel('SRT difference (predicted) [dB]')  
title('Prediction of SRT difference')
axis([-15 15 -15 15])
text(-14,6,['{\itN} = ' num2str(length(D_x))])

text(-14,14,['R^2 = ' num2str(round(r2,2))])
text(-14,12,['Bias = ' num2str(round(bias,2)) ' dB'])
text(-14,10,['RMSE = ' num2str(round(rmsv,2)) ' dB'])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figw figh]);

if sflag
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig5_scatter_SRTdiff_predicted.eps'],'-painters','-depsc','-r300');
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig5_scatter_SRTdiff_predicted.png'],'-dpng','-r300');
    
end
 

%% SRT difference (predicted) depending on (WRS_i - 50%), for different SII slopes 

beta0 = 1.522; 
beta1 = -0.515; 
beta2 = -0.146; 

s_ex = [1 2 3 4]; 
wrs_vec = 15:5:85; 

figh7 = figure('Visible','off'); 
hold on;
box on; 

line([15-50 85-50],[0 0],[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5])

for s = 1:length(s_ex)
    plot(wrs_vec-50,beta0+beta1*s_ex(s)+beta2*(wrs_vec-50),'linewidth',1.5) 
end

xlabel('WRS_i - 50 [%]')
ylabel('SRT difference (predicted) [dB]') 
legend('Perfect prediction','s_{SII} = 1 %/dB','s_{SII} = 2 %/dB','s_{SII} = 3 %/dB','s_{SII} = 4 %/dB') 
xlim([15 85]-50)
title(['GLM with \beta_0 = ' num2str(beta0) ', \beta_1 = ' num2str(beta1) ', \beta_2 = ' num2str(beta2) ] )



set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figw figh]);

if 0 %sflag
    print(figh7,[fig_folder filesep pfig_paper filesep 'Fig5B_SRTdiff_predicted_SIdiff.eps'],'-painters','-depsc','-r300');
    print(figh7,[fig_folder filesep pfig_paper filesep 'Fig5B_SRTdiff_predicted_SIdiff.png'],'-dpng','-r300');
    
end
