% Code for paper
% "Discrimination loss vs. SRT: A model-based approach towards harmonizing speech test interpretations"
% Buhl, Kludt, Schell-Majoor, Avan*, Campi*
% submitted to International Journal of Audiology
%
% Remote script 2: SII calculations
%
% Mareike Buhl
% mareike.buhl@pasteur.fr
%

clear all;
close all;
clc;

%% define output paths
dir_str = 'hl0-LD-mhh_opt_v2';
pfig = 'paper';

%% load individual audiograms

% change here to have all paths correct:
dataflag = 1; % 1: remote analysis (EK), 2: synthetic example analysis (MB)

warning('off');

if dataflag == 1
    addpath(genpath('./functions_remote/')); % all final functions need to be put here
    res_folder = './results_remote/';
    fig_folder = './figures_remote/';
elseif dataflag == 2
    path_str = '_remote'; % '_remote'|''
    addpath(genpath(['./functions' path_str '/']));
    res_folder = './results_example/';
    fig_folder = './figures_example/';
end

% load T_AG_SII from s1
load([res_folder pfig filesep 's1_T_AG_SII']);

ag_data_indiv = T_AG_SII{:,["id_patient","AC_250_Hz","AC_500_Hz","AC_1000_Hz","AC_1500_Hz","AC_2000_Hz","AC_3000_Hz","AC_4000_Hz","AC_6000_Hz","AC_8000_Hz"] };
ag_freq_indiv = [250 500 1000 1500 2000 3000 4000 6000 8000];

%% reference data
speechtest_refdata = [-7.1 19.9; ... % OLSA; SRT50 stationary noise and SRT50 in quiet (L0)
    -6.2 19.6; ... % GOESA
    -1.5 29.3]; % FMST

%% define conditions -> speechtest, ag data source, quiet/noise, ...
speechtests = {'FMST'};
quiet_noise = {'quiet'};
agtypes = {'_indiv'};

models = {'sii'};

calc_count = 0;
skip_count = 0;

for agtype = 1:length(agtypes)
    for m = 1:length(models)
        for cur_speechtest = 1:length(speechtests)
            for qn = 1:length(quiet_noise)
                disp([speechtests{cur_speechtest} ' - ' quiet_noise{qn} ' - ' agtypes{agtype} ' - ' models{m}]);

                %% from cfg script:
                targetdir = ['./' res_folder '/' pfig '/s2_' models{m} '/' dir_str];


                if ~exist([targetdir], 'dir')
                    mkdir([targetdir]);
                end

                % speechtest-specific:
                noisedir = '../../speechtests/noises';
                switch speechtests{cur_speechtest}
                    case 'FMST'
                        % reference parameters
                        speechref = {['../../speechtests/frei_add' filesep 'SAN-FBE.wav']};
                        noiseref = {[noisedir '/noises_ccitt.wav']};
                        genderref = 'male';

                        % speech and noise files
                        speeches = {['../../speechtests/frei_add' filesep 'SAN-FBE.wav']};
                        noises = {dir([noisedir filesep '*noises_ccitt.wav']).name,dir([noisedir filesep '*icra5_0-25_fullscale.wav']).name};
                        gender = {'male'};
                end


                switch quiet_noise{qn}
                    case 'noise'
                        srtref = speechtest_refdata(cur_speechtest,1);
                        switch cur_speechtest
                            case 3
                                lnlevelref = 60;
                            otherwise
                                lnlevelref = 65;
                        end
                        lslevelref = lnlevelref+srtref;
                        SNR = -25:5:60; % complete SNR vector - -25 instead of -15 für ICRA5-250
                        noise_levels = lnlevelref;

                    case 'quiet'
                        srtref = speechtest_refdata(cur_speechtest,2);
                        lnlevelref = -50; % -> -50 dB according to the standard
                        lslevelref = srtref;
                        SNR = 40:5:170; % complete SNR vector
                        noise_levels = -50;

                        noises = noises(1); % not used for SII in quiet (dummy variable for the loop)
                end

                language = {'DE'};

                % Audiogram
                if strcmp(agtypes{agtype},'_bisg')
                    load('./results/data_lit/bisgaard_ags.mat'); % variables bis_ags, bis_names, bis_freqs
                    bis_ags =  [0 0 0 0 0 0 0 0 0 0; ... % NH
                        bis_ags];
                    hloss_ID = [{'NH'},bis_names];

                    hloss = bis_ags;

                    % calculate reference SII
                    disp('ref-sii')
                    FNH = bis_freq;
                    NH = hloss(1,:);
                elseif strcmp(agtypes{agtype},'_indiv') % use loaded data for individual patients
                    hloss = ag_data_indiv(:,2:2+length(ag_freq_indiv)-1); % first column contains id (MHH)
                    hloss_ID = ag_data_indiv(:,1);
                    FNH = ag_freq_indiv;
                    NH = zeros(size(ag_freq_indiv));
                end


                LSref = srtref;
                LNref = lnlevelref;
                switch quiet_noise{qn}
                    case 'noise'
                        LSref = srtref + LNref; % for SII in noise, distinction only needed for reference value (because to-be-included SRT/speech level is given)
                    case 'quiet'
                        LSref = srtref; % in quiet
                end

                tmp = strsplit(noiseref{1},'/');
                tmp2 = strsplit(tmp{end},'.');


                noiseref_name = tmp2{1};
                if strcmp(noiseref_name,'noises_ccitt')
                    noiseref_name = 'ccitt';
                elseif strcmp(noiseref_name,'icra5_0-25_fullscale')
                    noiseref_name = 'icra5250';
                end

                cmv(models{m},speechref,noiseref{1},srtref,genderref,[],[],NH,FNH,[targetdir filesep models{m} '_' speechtests{cur_speechtest} '_' strrep(noiseref_name,'-','_') '_' num2str(LNref) '_' quiet_noise{qn} '-ref.txt'],LSref,LNref); % 02/01/23


                %% calculate desired model values
                T_slope = table();
                fname_flag = 0;

                for lnlevel = noise_levels
                    for ilanguage = 1
                        for knoise = 1:length(noises)
                            tmp = strsplit(noises{knoise},'.');
                            noise_name = tmp{1};
                            if strcmp(noise_name,'noises_ccitt')
                                noise_name = 'ccitt';
                            elseif strcmp(noise_name,'icra5_0-25_fullscale')
                                noise_name = 'icra5250';
                            end
                            sfile = speeches(ilanguage);
                            for jhl = 1:length(hloss_ID)
                                if mod(jhl,100) == 0
                                    disp(['Progress: ' num2str(round(jhl/length(hloss_ID)*100,1)) ' %']);
                                end
                                if isnumeric(hloss_ID)
                                    hlIDstr = ['p' num2str(hloss_ID(jhl))];
                                    T_slope.id_patient(jhl) = {hlIDstr};
                                else
                                    hlIDstr = hloss_ID{jhl};
                                    T_slope.id_patient(jhl) = {hlIDstr};
                                end

                                level_vec = lnlevel+SNR;
                                sii_vec = nan(size(level_vec));

                                level_start = [5 11 17 23];
                                levels = level_start;

                                for lev_idx = 1:length(levels)
                                    lslevel = level_vec(levels(lev_idx));

                                    fname = [targetdir filesep models{m} '_' speechtests{cur_speechtest} '_' hlIDstr '_' language{ilanguage} '_' num2str(lslevel) '_' strrep(noise_name,'-','_') '_' num2str(lnlevel) '_' quiet_noise{qn} agtypes{agtype} '.txt'];

                                    if  ~isfile(fname) % check if condition was already calculated
                                        sii_vec(levels(lev_idx)) = cmv(models{m},sfile,[noisedir filesep noises{knoise}],srtref,gender{ilanguage},[],[],hloss(jhl,:),FNH,fname,lslevel,lnlevel);
                                        calc_count = calc_count + 1;
                                        sii_vec(levels(lev_idx));
                                        fname_flag = 1;
                                    else
                                        skip_count = skip_count + 1;
                                    end


                                end % speech level

                                if  fname_flag
                                    fname_flag = 0; % to conduct this part always if SII calculation is conducted above

                                    % estimate additional sii values until a reliable slope is
                                    % obtained
                                    r2_max = 0; % to get started
                                    step = 1;
                                    while r2_max < 0.99
                                        % estimate additional levels
                                        [levels,level_add] = add_levels(levels,sii_vec,step);
                                        if isempty(level_add)
                                            break;
                                        end

                                        % calculate SII for new levels
                                        for lev_idx = 1:length(level_add)
                                            sii_vec(level_add(lev_idx)) = cmv(models{m},sfile,[noisedir filesep noises{knoise}],srtref,gender{ilanguage},[],[],hloss(jhl,:),FNH,fname,level_vec(level_add(lev_idx)),lnlevel);
                                        end

                                        % estimate if obtained linear area can be used for slope
                                        % estimation
                                        sii_tmp = sii_vec(levels);
                                        [r2_max,r2_idx] = check_lin_slope(levels,sii_tmp);
                                        step = step+1;
                                    end
                                    if isempty(level_add)
                                        T_slope.slope(jhl) = nan;
                                        T_slope.r2(jhl) = nan;
                                    else
                                        T_slope.slope(jhl) = (sii_tmp(r2_idx+1)-sii_tmp(r2_idx-1))/(level_vec(levels(r2_idx+1))-level_vec(levels(r2_idx-1))); % would for the moment not save if more than one noise condition or language
                                        T_slope.r2(jhl) = r2_max;
                                    end
                                    clear sii_tmp;

                                    % save slope results
                                    writetable(T_slope(jhl,:),[targetdir filesep 'Tslope_' models{m} '_' speechtests{cur_speechtest} '_' language{ilanguage} '_' strrep(noise_name,'-','_') '_' num2str(lnlevel) '_' quiet_noise{qn} agtypes{agtype} '_' hlIDstr],FileType='text')
                                end

                            end
                        end
                    end
                end
            end
        end
    end
end

