% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Local script 6: Plot figure 6
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

filter_crit_red_all = (T_fit_emp_all_sii.case == 3 | T_fit_emp_all_sii.case == 5  ) & T_fit_emp_all_sii.srt_pta_diff_80 >= thr_c;  

%% Figure 6 

D_x1 = 100-T_fit_emp_all_sii.maxSI(filter_crit_red_all);
D_y1 = T_fit_emp_all_sii.D_srt_f1(filter_crit_red_all);   

% markersize depending on frequency of same data point 
[D_unique,ia,ic] = unique([5*round(D_x1/5),2*round(D_y1/2)],'rows'); 
h_uni1 = hist(ic,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni1 = h_uni1/sum(h_uni1); 

[P,xbins] = calc_percentiles_sections(100-T_fit_emp_all_sii.maxSI,T_fit_emp_all_sii.D_srt_f1,0:5:40,filter_crit_red_all,10:20:90); %0:10:100


D_x2 = 100-T_fit_emp_all_sii.maxSI(filter_crit_yellow_all); % filter_crit_yellow_all
D_y2 = T_fit_emp_all_sii.D_srt_fsii(filter_crit_yellow_all);   

% markersize depending on frequency of same data point 
[D_unique2,ia,ic2] = unique([5*round(D_x2/5),2*round(D_y2/2)],'rows'); 
h_uni2 = hist(ic2,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni2 = h_uni2/sum(h_uni2); 

[P2,xbins] = calc_percentiles_sections(100-T_fit_emp_all_sii.maxSI,T_fit_emp_all_sii.D_srt_fsii,0:5:40,filter_crit_yellow_all,10:20:90); 

D_x3 = 100-T_fit_emp_all_sii.maxSI(filter_crit_violet_all);
D_y3 = T_fit_emp_all_sii.D_srt_fall(filter_crit_violet_all);   

% markersize depending on frequency of same data point 
[D_unique3,ia,ic3] = unique([5*round(D_x3/5),2*round(D_y3/2)],'rows'); 
h_uni3 = hist(ic3,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni3 = h_uni3/sum(h_uni3); 

[P3,xbins] = calc_percentiles_sections(100-T_fit_emp_all_sii.maxSI,T_fit_emp_all_sii.D_srt_fall,0:5:40,filter_crit_violet_all,10:20:90); 

figh2 = figure; 

% violet
subplot(1,3,1)
hold on;
box on;

% "legend"
legend_size3 = [round(max(h_uni3)*length(ic3)) round(max(h_uni3)*length(ic3)/2) round(max(h_uni3)*length(ic3)/4) round(max(h_uni3)*length(ic3)/8) round(max(h_uni3)*length(ic3)/16)]; 
legend_y3 = 47; 
for n = 1:length(legend_size3)
    scatter(32,legend_y3-(n-1)*3,4000*log10(legend_size3(n)/length(ic3)+1),'MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:));  
    text(34,legend_y3-(n-1)*3,num2str(round(legend_size3(n))))
end 

scatter(D_unique3(:,1),D_unique3(:,2),4000*log10(h_uni3+1)','MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:)); %,'o','linewidth',1.5,'color',fig_properties.fillcolors(3,:),'MarkerFaceColor',fig_properties.fillcolors(3,:));

% percentiles 
lwidth = [0.5 1.5 2.5 1.5 0.5  ]; 
for n = 1:size(P3,2)
    plot(xbins,P3(:,n),'-','color',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.fillcolors(4,:),'linewidth',lwidth(n)); 
end

ylabel('D_{SRT,f} [dB]')
xlabel('100 - WRS_{max} [%]')
title('Empirical slope')
axis([-2 42 -10 50])
text(1,48,['{\itN} = ' num2str(length(ic3))])


% yellow
subplot(1,3,2)
hold on;
box on;

% "legend"
legend_size2 = [round(max(h_uni2)*length(ic2)) round(max(h_uni2)*length(ic2)/2) round(max(h_uni2)*length(ic2)/4) round(max(h_uni2)*length(ic2)/8) round(max(h_uni2)*length(ic2)/16)]; 
legend_y2 = 47; 
for n = 1:length(legend_size2)
    scatter(32,legend_y2-(n-1)*3,4000*log10(legend_size2(n)/length(ic2)+1),'MarkerEdgeColor',fig_properties.fillcolors(3,:),'MarkerFaceColor',fig_properties.linecolors2(3,:));  
    text(34,legend_y2-(n-1)*3,num2str(round(legend_size2(n))))
end 

scatter(D_unique2(:,1),D_unique2(:,2),4000*log10(h_uni2+1)','MarkerEdgeColor',fig_properties.fillcolors(3,:),'MarkerFaceColor',fig_properties.linecolors2(3,:));

% percentiles 
for n = 1:size(P2,2)
    plot(xbins,P2(:,n),'-','color',fig_properties.fillcolors(3,:),'MarkerFaceColor',fig_properties.fillcolors(3,:),'linewidth',lwidth(n)); 
end

ylabel('D_{SRT,h} [dB]')
xlabel('100 - WRS_{max} [%]')
title('SII-based slope')
axis([-2 42 -10 50])
text(1,48,['{\itN} = ' num2str(length(ic2))])


% red 
subplot(1,3,3)
hold on;
box on;

% "legend"
legend_size1 = [round(max(h_uni1)*length(ic)) round(max(h_uni1)*length(ic)/2) round(max(h_uni1)*length(ic)/4) round(max(h_uni1)*length(ic)/8) round(max(h_uni1)*length(ic)/16)]; 
legend_y1 = 46; 
for n = 1:length(legend_size1)
    scatter(30,legend_y1-(n-1)*5,4000*log10(legend_size1(n)/length(ic)+1),'MarkerEdgeColor',fig_properties.fillcolors(2,:),'MarkerFaceColor',fig_properties.linecolors2(2,:));  
    text(34,legend_y1-(n-1)*5,num2str(round(legend_size1(n))))
end 

scatter(D_unique(:,1),D_unique(:,2),4000*log10(h_uni1+1)','MarkerEdgeColor',fig_properties.fillcolors(2,:),'MarkerFaceColor',fig_properties.linecolors2(2,:)); 

% percentiles 
for n = 1:size(P,2)
    plot(xbins,P(:,n),'-','color',fig_properties.fillcolors(2,:),'MarkerFaceColor',fig_properties.fillcolors(2,:),'linewidth',lwidth(n)); 
end

ylabel('D_{SRT,n} [dB]')
xlabel('100 - WRS_{max} [%]')
title('Normal-hearing slope')
axis([-2 42 -10 50])
text(1,48,['{\itN} = ' num2str(length(ic))])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.5*figw figh]);

if sflag
    print(figh2,[fig_folder filesep pfig_paper filesep 'Fig6_Dsrt_maxSI.eps'],'-painters','-depsc','-r300');
    print(figh2,[fig_folder filesep pfig_paper filesep 'Fig6_Dsrt_maxSI.png'],'-dpng','-r300');
    
end

