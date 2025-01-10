
% INPUT
% num_params    desired number and choice of fit parameters (1: SRT; 2: SRT
%               and slope; 3: SRT, slope, maxSI)

function [srt,slp,maxsi,gof] = fit_srt_psyfun(level_vec,si_vec,num_params)

srt = nan;
slp = nan;
maxsi = nan;
gof = struct('sse', {nan}, 'rsquare', {nan}, 'dfe', {nan},'adjrsquare', {nan},'rmse', {nan});

if ~isnan(level_vec)

    if num_params == 1
        if length(si_vec) >= 1
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[25],... % 0
                'Upper',[150],...
                'StartPoint',[70]);

            if si_vec == 100
                level_vec = [0;level_vec]; % needed for best performing patients (SI=100%)
                si_vec = [0;si_vec];
            end

            [f,gof,output] = fit(level_vec,si_vec,'100/(1+exp(4*4.5/100*(b-x)))',fo);
            % requires column vectors
            srt = f.b;
            slp = nan;
            maxsi = nan;
        end

    elseif num_params == 2
        if length(si_vec) >= 2
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0, 20],...
                'Upper',[0.05,150],... % slope max different to v1 of fit
                'StartPoint',[0.05, 70]);

            if si_vec == 100
                level_vec = [0;level_vec];
                si_vec = [0;si_vec];
            end

            [f,gof,output] = fit(level_vec,si_vec,'100/(1+exp(4*a*(b-x)))',fo);
            srt = f.b;
            slp = f.a;
            maxsi = nan;
        end

    elseif num_params == 3
        if length(si_vec) >= 3
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0, 0, 0],...
                'Upper',[0.05,150,100],... % 0.05
                'StartPoint',[0.05, 70, 100]);

            [f,gof,output] = fit(level_vec,si_vec,'c/(1+exp(4*a*(b-x)))',fo);
            srt = f.b;
            slp = f.a;
            maxsi = f.c;

        end

    end
end

