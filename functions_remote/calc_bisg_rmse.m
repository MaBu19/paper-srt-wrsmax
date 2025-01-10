% calc_bisg_rmse (RMS error between two audiograms -> e.g. one individual and
% one Bisgaard audiogram)
%
% Mareike Buhl
% mareike.buhl@pasteur.fr
%
% v1.0, 10.08.2022
% v2.0, 16.07.2024: added class 'none' if no input data 
%
% INPUT
% ag_data             vector containing row of audiograms 
% data_freq           frequency vector corresponding to ag_data 
%
% OUTPUT
% bis_class           label of respective Bisgaard class (important when
%                     working with tables, e.g. splitapply())
% ag_rmse_out         RMS error of closest Bisgaard audiogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bis_class ag_rmse_out] = calc_bisg_rmse(ag_data, data_freq)

% load Bisgaard audiograms
load('./results/data_lit/bisgaard_ags.mat');
bis_names = categorical(bis_names);

% estimate common frequencies (could also be interpolated)
[common_freq,idx_data_freq,idx_bis_freq] = intersect(data_freq,bis_freq);

% calculate RMS error (difference)
idx_nan = ~isnan(ag_data(:,idx_data_freq));
ag_dist = bis_ags(:,idx_bis_freq(idx_nan))-ag_data(:,idx_data_freq(idx_nan));
ag_rmse = rms(ag_dist',1);
if ~isnan(ag_rmse)
    [~,id_bis_next] = min(ag_rmse);
    ag_rmse_out = ag_rmse(id_bis_next);
    bis_class = bis_names(id_bis_next);
else
    ag_rmse_out = nan;
    bis_class = 'none';
end
end