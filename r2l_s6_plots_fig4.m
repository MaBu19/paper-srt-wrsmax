% Code for paper 
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations" 
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology 
% 
% Local script 6: Plot figure 4  
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
subflag = 1; % 1 subplot, 2 separate plots 


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

%% Figure 4: Scatterplots SRT, slope, D

% % calculate percentiles: 
[P,xbins] = calc_percentiles_sections(T_fit_emp_all_sii.fall_srt-srt_ref,T_fit_emp_all_sii.fsii_indiv_srt-srt_ref,20:2.5:60,filter_crit_violet,10:20:90);

% Figure 3 A: SRT_SII vs SRT_DD 
D_x = T_fit_emp_all_sii.fall_srt(filter_crit_violet) - srt_ref; 
D_y = T_fit_emp_all_sii.fsii_indiv_srt(filter_crit_violet) - srt_ref; 

% markersize depending on frequency of same data point 
[D_unique,ia,ic] = unique([2*round(D_x/2),2*round(D_y/2)],'rows'); 
h_uni1 = hist(ic,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni1 = h_uni1/sum(h_uni1); 

plot_colors = 4*ones(size(D_x));  

idx_D = 1:length(D_x); 

[r2,p,bias,rmsv,rl2,ru2] = statistical_analysis(D_x(idx_D),D_y(idx_D));

if subflag == 2
    figh6 = figure('Visible','on');
elseif subflag == 1
    figh7 = figure('Visible','on');
    subplot(1,3,1)
end
hold on;
box on;

line([-5 120],[-5 120],[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

for i = 1:size(D_x)
    p(i) = plot(150,150,'o','linewidth',1.5,'color',fig_properties.linecolors2(plot_colors(i),:),'MarkerFaceColor',fig_properties.linecolors2(plot_colors(i),:)); % dummy plot for legend
    pleg(plot_colors(i)) = p(i);  
end

% "legend"
legend_size1 = [round(max(h_uni1)*length(ic)) round(max(h_uni1)*length(ic)/2) round(max(h_uni1)*length(ic)/4) round(max(h_uni1)*length(ic)/8) round(max(h_uni1)*length(ic)/16)]; 
legend_y1 = 32; 
for n = 1:length(legend_size1)
    scatter(52,legend_y1-(n-1)*3,4000*log10(legend_size1(n)/length(ic)+1),'MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:));  
    text(55,legend_y1-(n-1)*3,num2str(round(legend_size1(n))))
end 

lwidth = [0.5 1.5 2.5 1.5 0.5  ];
for n = 1:size(P,2)
    plot(xbins,P(:,n),'-','color',fig_properties.fillcolors(plot_colors(i),:),'MarkerFaceColor',fig_properties.fillcolors(plot_colors(i),:),'linewidth',lwidth(n))

end

scatter(D_unique(:,1),D_unique(:,2),4000*log10(h_uni1+1)','MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:)); 

xlabel('SRT_f - SRT_{NH} [dB]')
ylabel('SRT_h - SRT_{NH} [dB]')  
title('SRT loss')
axis([20 58 20 58])
text(22,49,['{\itN} = ' num2str(length(D_x))])

text(22,56,['R^2 = ' num2str(round(r2,2))])
text(22,54,['Bias = ' num2str(round(bias,2)) ' dB'])
text(22,52,['RMSE = ' num2str(round(rmsv,2)) ' dB'])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figw figh]);

if subflag == 2
if sflag
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig4A_scatter_SRTloss_twomodes.eps'],'-painters','-depsc','-r300');
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig4A_scatter_SRTloss_twomodes.png'],'-dpng','-r300');
    
end
end 

%% Figure 3 B: slope_SII vs slope_DD 

[Ps,xbins] = calc_percentiles_sections(T_fit_emp_all_sii.fall_slope,T_fit_emp_all_sii.fsii_indiv_slope,0:0.25:4.5,filter_crit_violet,10:20:90); 
 
D_x = T_fit_emp_all_sii.fall_slope(filter_crit_violet); 
D_y = T_fit_emp_all_sii.fsii_indiv_slope(filter_crit_violet);  
plot_colors = 4*ones(size(D_x));    

% markersize depending on frequency of same data point 
[D_unique,ia,ic] = unique([0.25*round(D_x/0.25),0.25*round(D_y/0.25)],'rows'); 
h_uni1 = hist(ic,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni1 = h_uni1/sum(h_uni1); 

idx_D = 1:length(D_x);

[r2,p,bias,rmsv,rl2,ru2] = statistical_analysis(D_x,D_y);

if subflag == 2
    figh6 = figure('Visible','on');
elseif subflag == 1
    subplot(1,3,2)
end
hold on;
box on;

line([0 4.6],[0 4.6],[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

% "legend"
legend_size1 = [round(max(h_uni1)*length(ic)) round(max(h_uni1)*length(ic)/2) round(max(h_uni1)*length(ic)/4) round(max(h_uni1)*length(ic)/8) round(max(h_uni1)*length(ic)/16)]; 
legend_y1 = 3; 
for n = 1:length(legend_size1)
    scatter(3.5,legend_y1-(n-1)*0.3,4000*log10(legend_size1(n)/length(ic)+1),'MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:));  
    text(3.7,legend_y1-(n-1)*0.3,num2str(round(legend_size1(n))))
end 

lwidth = [0.5 1.5 2.5 1.5 0.5  ];
for n = 1:size(Ps,2)
    plot(xbins,Ps(:,n),'-','color',fig_properties.fillcolors(plot_colors(i),:),'MarkerFaceColor',fig_properties.fillcolors(plot_colors(i),:),'linewidth',lwidth(n))

end

scatter(D_unique(:,1),D_unique(:,2),4000*log10(h_uni1+1)','MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:)); 

xlabel('Slope s_{f} [%/dB]')
ylabel('Slope s_{h} [%/dB]') 
title('Slope')
axis([0 4.6 0 4.6])

text(2.8,1,['R^2 = ' num2str(round(r2,2))])
text(2.8,0.7,['Bias = ' num2str(round(bias,2)) ' dB'])
text(2.8,0.4,['RMSE = ' num2str(round(rmsv,2)) ' dB'])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figw figh]);

if subflag == 2
if sflag
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig4B_scatter_slope_twomodes.eps'],'-painters','-depsc','-r300');
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig4B_scatter_slope_twomodes.png'],'-dpng','-r300');
    
end
end 


%% Figure 3 C: D_SII vs D_DD 

[PD,xbins] = calc_percentiles_sections(T_fit_emp_all_sii.D_srt_fall,T_fit_emp_all_sii.D_srt_fsii,-10:2.5:50,filter_crit_violet,10:20:90); 

D_x = T_fit_emp_all_sii.D_srt_fall(filter_crit_violet); 
D_y = T_fit_emp_all_sii.D_srt_fsii(filter_crit_violet); 

% markersize depending on frequency of same data point 
[D_unique,ia,ic] = unique([2*round(D_x/2),2*round(D_y/2)],'rows'); 
h_uni1 = hist(ic,1:length(ia)); % ia contains corresponding data points in the same order 
h_uni1 = h_uni1/sum(h_uni1); 

plot_colors = 4*ones(size(D_x)); 
idx_D = 1:length(D_x);

[r2,p,bias,rmsv,rl2,ru2] = statistical_analysis(D_x(idx_D),D_y(idx_D));

if subflag == 2
    figh6 = figure('Visible','on');
elseif subflag == 1
    subplot(1,3,3)
end
hold on;
box on;

line([-10 120],[-10 120],[0 0],'linewidth',1.5,'color',1.5*[0.5 0.5 0.5]) 

% "legend"
legend_size1 = [round(max(h_uni1)*length(ic)) round(max(h_uni1)*length(ic)/2) round(max(h_uni1)*length(ic)/4) round(max(h_uni1)*length(ic)/8) round(max(h_uni1)*length(ic)/16)]; 
legend_y1 = 8; 
for n = 1:length(legend_size1)
    scatter(40,legend_y1-(n-1)*3,4000*log10(legend_size1(n)/length(ic)+1),'MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:));  
    text(42,legend_y1-(n-1)*3,num2str(round(legend_size1(n))))
end 

lwidth = [0.5 1.5 2.5 1.5 0.5  ];
for n = 1:size(PD,2)
    plot(xbins,PD(:,n),'-','color',fig_properties.fillcolors(plot_colors(i),:),'MarkerFaceColor',fig_properties.fillcolors(plot_colors(i),:),'linewidth',lwidth(n))

end

scatter(D_unique(:,1),D_unique(:,2),4000*log10(h_uni1+1)','MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:)); 


xlabel('D_{SRT,f} [dB]') 
ylabel('D_{SRT,h} [dB]') 
title('D component')
axis([-10 50 -10 50])

text(-7,47,['R^2 = ' num2str(round(r2,2))])
text(-7,43,['Bias = ' num2str(round(bias,2)) ' dB'])
text(-7,39,['RMSE = ' num2str(round(rmsv,2)) ' dB'])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figw figh]);

if subflag == 2
if sflag
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig4C_scatter_Dsrt_twomodes.eps'],'-painters','-depsc','-r300');
    print(figh6,[fig_folder filesep pfig_paper filesep 'Fig4C_scatter_Dsrt_twomodes.png'],'-dpng','-r300');
    
end
end 

if subflag == 1
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.5*figw figh]);

if sflag
    print(figh7,[fig_folder filesep pfig_paper filesep 'Fig4_comparison_subplot.eps'],'-painters','-depsc','-r300');
    print(figh7,[fig_folder filesep pfig_paper filesep 'Fig4_comparison_subplot.png'],'-dpng','-r300');
    
end
end 
