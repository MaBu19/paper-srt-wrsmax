function [itc,slp,srt,gof] = fit_linear_slope(level_vec,si_vec,num_params,sii_slope)

itc = nan;
slp = nan;
srt = nan;
gof = struct('sse', {nan}, 'rsquare', {nan}, 'dfe', {nan},'adjrsquare', {nan},'rmse', {nan});

xvec = 10:110; % for SRT estimation
 
if num_params == 1 % linear fit based on min. 2 data points
    if length(si_vec) >= 2 && any(diff(si_vec)>0)
        fo = fitoptions('Method','NonlinearLeastSquares');

        if si_vec == 100
            level_vec = [0;level_vec]; % needed for best performing patients (SI=100%)
            si_vec = [0;si_vec];
        end

        [f,gof,output] = fit(level_vec,si_vec,'a*x+b',fo);
        % requires column vectors
        itc = f.b;
        slp = f.a;

        % SRT from linear function:
        srt = interp1(itc+slp*xvec,xvec,50);

    end

elseif num_params == 2 % linear fit based on 1 data point and given slope
    if ~isnan(sii_slope)

        if length(si_vec) >= 1
            fo = fitoptions('Method','NonlinearLeastSquares');

            if si_vec == 100
                level_vec = [0;level_vec]; % needed for best performing patients (SI=100%)
                si_vec = [0;si_vec];
            end

            eq = 'a*x+b';
            % include fixed slope as given by SII (relative to NH slope)
            eq = strrep(eq,'a',num2str(sii_slope));

            [f,gof,output] = fit(level_vec,si_vec,eq,fo);
            % requires column vectors
            itc = f.b;
            slp = sii_slope; % same as input

            % SRT from linear function:
            if slp > 0 % otherwise sample points not unique
                srt = interp1(itc+slp*xvec,xvec,50);
            end

        end
    end
end
end


