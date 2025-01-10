% Code for paper
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations"
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology
%
% Local script 6: Plot figure 2
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
fig_properties.fillcolors(2,:) = [255 114 118]/255;
fig_properties.fillcolors(3,:) = [241 235 156]/255;
fig_properties.fillcolors(4,:) = [177 162 250]/255;


%% filter criteria
filter_crit_violet = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2 & ~isnan(T_fit_emp_all_sii.srt_diff) & T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI;
filter_crit_violet_all = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 2 & T_fit_emp_all_sii.srt_pta_diff_dd >= thr_c;

filter_crit_yellow_all = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI;
filter_crit_yellow_1DP = T_fit_emp_all_sii.case == 2 & T_fit_emp_all_sii.m_dp_fit == 1 & T_fit_emp_all_sii.srt_pta_diff_sii >= thr_c & T_fit_emp_all_sii.SI_sii >= T_fit_emp_all_sii.maxSI;

filter_crit_red_all = (T_fit_emp_all_sii.case == 3  | T_fit_emp_all_sii.case == 5) & T_fit_emp_all_sii.srt_pta_diff_80 >= thr_c;

filter_crit_blue_all = T_fit_emp_all_sii.case == 4 | T_fit_emp_all_sii.case == 6;

%% Figure A: stacked histogram of available measured levels

if subflag == 2
    figh5 = figure('Visible','on');
elseif subflag == 1
    figh7 = figure('Visible','on');
    subplot(1,2,1)
end

hold on; box on;

h1 = hist(T_fit_emp_all_sii.m_dp(filter_crit_blue_all),[0 1 2 3]);  
h2 = hist(T_fit_emp_all_sii.m_dp(filter_crit_red_all),[0 1 2 3]);  
h3 = hist(T_fit_emp_all_sii.m_dp(filter_crit_yellow_1DP),[0 1 2 3]); 
h4 = hist(T_fit_emp_all_sii.m_dp(filter_crit_violet_all),[0 1 2 3]); 

b1 = bar([1 2 3],[h1(2:4); h2(2:4); h3(2:4); h4(2:4)]','stacked');

xlabel('Number of measured test lists')
ylabel('Number of patients')
title('Data availability for SRT estimation')
axis([0.5 3.5 0 max(sum([h1; h2; h3; h4]))*1.1])
text(2.3,0.68*max(sum([h1; h2; h3; h4]))*0.97,['{\itN} = ' num2str(sum(sum([h1; h2; h3; h4])))])
legend('No SRT feasible','Undetermined','Half-determined','Fully-determined','Location','NorthEast')
set(gca,'XTick',1:1:3)

if subflag == 2
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figw figh]);

    if sflag
        print(figh5,[fig_folder filesep pfig_paper filesep 'Fig1A_histo_DP_cases.eps'],'-painters','-depsc','-r300');
        print(figh5,[fig_folder filesep pfig_paper filesep 'Fig1A_histo_DP_cases.png'],'-dpng','-r300');

    end
end

%% Figure B: maxSI vs. PTA

D_x = T_fit_emp_all_sii.pta_spl; 
D_y = T_fit_emp_all_sii.maxSI;
idx_D = ones(size(D_x));
scaling = 30; % for adjusting marker size in plot - normalization needed for relative size between cases

[D_unique1,ia1,ic1] = unique([5*round(T_fit_emp_all_sii.pta_spl(filter_crit_blue_all)/5),5*round(T_fit_emp_all_sii.maxSI(filter_crit_blue_all)/5)],'rows');
h_uni1 = hist(ic1,1:length(ia1)); % ia contains corresponding data points in the same order
h_uni1 = h_uni1/sum(h_uni1); 

[D_unique2,ia2,ic2] = unique([5*round(T_fit_emp_all_sii.pta_spl(filter_crit_red_all)/5),5*round(T_fit_emp_all_sii.maxSI(filter_crit_red_all)/5)],'rows');
h_uni2 = hist(ic2,1:length(ia2)); 
h_uni2 = h_uni2/sum(h_uni2);

[D_unique3,ia3,ic3] = unique([5*round(T_fit_emp_all_sii.pta_spl(filter_crit_yellow_1DP)/5),5*round(T_fit_emp_all_sii.maxSI(filter_crit_yellow_1DP)/5)],'rows');
h_uni3 = hist(ic3,1:length(ia3)); 
h_uni3 = h_uni3/sum(h_uni3);

[D_unique4,ia4,ic4] = unique([5*round(T_fit_emp_all_sii.pta_spl(filter_crit_violet_all)/5),5*round(T_fit_emp_all_sii.maxSI(filter_crit_violet_all)/5)],'rows');
h_uni4 = hist(ic4,1:length(ia4)); 
h_uni4 = h_uni4/sum(h_uni4);

if subflag == 2
    figh6 = figure('Visible','on');
elseif subflag == 1
    subplot(1,2,2)
end
hold on;
box on;

p = scatter([150 150 150]',[150 150 150]',70*log10([1 10 100]'+1),'k');

% "legend"
legend_size1 = [round(max(h_uni4)*length(ic4)) round(max(h_uni4)*length(ic4)/2) round(max(h_uni4)*length(ic4)/4) round(max(h_uni4)*length(ic4)/8) round(max(h_uni4)*length(ic4)/16)];
legend_y1 = 30;
for n = 1:length(legend_size1)
    scatter(0,legend_y1-(n-1)*5,scaling*70*log10(legend_size1(n)/length(ic1)+1),'MarkerEdgeColor',fig_properties.fillcolors(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:));
    text(4,legend_y1-(n-1)*5,num2str(round(legend_size1(n))))
end

scatter(D_unique1(:,1)-1.5,D_unique1(:,2),scaling*70*log10(h_uni1+1)','MarkerEdgeColor',fig_properties.linecolors2(1,:),'MarkerFaceColor',fig_properties.linecolors2(1,:));
scatter(D_unique2(:,1)-0.5,D_unique2(:,2),scaling*70*log10(h_uni2+1)','MarkerEdgeColor',fig_properties.linecolors2(2,:),'MarkerFaceColor',fig_properties.linecolors2(2,:));
scatter(D_unique3(:,1)+0.5,D_unique3(:,2),scaling*70*log10(h_uni3+1)','MarkerEdgeColor',fig_properties.linecolors2(3,:),'MarkerFaceColor',fig_properties.linecolors2(3,:));
scatter(D_unique4(:,1)+1.5,D_unique4(:,2),scaling*70*log10(h_uni4+1)','MarkerEdgeColor',fig_properties.linecolors2(4,:),'MarkerFaceColor',fig_properties.linecolors2(4,:));


xlabel(['PTA [dB SPL]'])
ylabel(['WRS_{max} [%]'])
axis([-5 120 0 103])
title('Empirical data')

if subflag == 2
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figw figh]);

    if sflag
        print(figh6,[fig_folder filesep pfig_paper filesep 'Fig1B_scatter_maxSI_PTA_colorcases.eps'],'-painters','-depsc','-r300');
        print(figh6,[fig_folder filesep pfig_paper filesep 'Fig1B_scatter_maxSI_PTA_colorcases.png'],'-dpng','-r300');

    end
end

if subflag == 1
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 2.3*figw figh]);

    if sflag
        print(figh7,[fig_folder filesep pfig_paper filesep 'Fig2_data_subplot.eps'],'-painters','-depsc','-r300');
        print(figh7,[fig_folder filesep pfig_paper filesep 'Fig2_data_subplot.png'],'-dpng','-r300');

    end
end

