function mv = cmv(model,speech_file,noise_file,SRTref,gender,mvref,dmvref,hlloss,hlfreq,target_file,LS,LN)

	bCalcSRT = true;
	if isempty(mvref)
		bCalcSRT = false;
	end

	% Parameters
	fs 			   = 44100;
	iterations = 1;       % run every model "iterations" times - 16
	speech_dur = 10; 			 %10 s
	fade_dur 	 = 10*10^-3; %10 ms

	[n,nfs] = audioread(noise_file);
	if max(size(speech_file)) == 1
		[s,sfs] = audioread(speech_file{1});

		if length(s) < sfs*speech_dur
			s = crossfade_extend(s, ceil(fade_dur*sfs), sfs*speech_dur);
		end

	else % e.g. used for STOI
		[s,sfs] = loadsentences(speech_file{1},speech_file{2});
	end
	if sfs ~= fs; s = resample(s, fs, sfs); end
	if nfs ~= fs; n = resample(n, fs, nfs); end 

    if min(size(s))
        s = s(:,1);
    end

	% prepare data
if length(SRTref) == 1 

[s, n] = adjSpeechNoiseLevel_refNoise(s, n, LS-LN, LN); 

elseif length(SRTref) == 2 
  [s, n] = adjSpeechNoiseLevel_refNoise(s, n, SRTref(1)-SRTref(2),SRTref(2));
  SRTref(1)-SRTref(2)
end


	% just calculate reference model value and return
	if ~bCalcSRT
        parfor step = 1:iterations 
            if length(s) > sfs*speech_dur 
                start_idx = floor(rand(1)*(length(s)-1.5*sfs*speech_dur));
                start_idx = max(0,start_idx); 
                s_rd = s(start_idx+1:start_idx+(sfs*speech_dur),:);

            else 
                s_rd = s; 
            end

			[signal, noise] = res(s_rd,n); 
			mvout_tmp(step) = calculate_objective_measures(signal,noise,fs,gender,hlloss,hlfreq,model);
        end

		mv   = mean(mvout_tmp);
		dmv  = std(mvout_tmp);
		SRT  = nan;
		dSRT = nan;
	end

	% just calculate SRT 
	if bCalcSRT
		tmpSRT = nan(iterations,1); tmpmv = nan(iterations,1);
		parfor jj = 1:iterations
			[signal, noise] = res(s,n);
      if length(SRTref) == 1 
			[tmpSRT(jj) tmpmv(jj)] = SRT_from_model_fast(signal, noise, fs, model, SRTref, mvref ,max(dmvref,0.001), hlloss, hlfreq, gender);
      elseif length(SRTref) == 2 
      			[tmpSRT(jj) tmpmv(jj)] = SRT_from_model_fast(signal, noise, fs, model, SRTref(1)-SRTref(2), mvref ,max(dmvref,0.001), hlloss, hlfreq, gender);
      end

		end
		SRT  = mean(tmpSRT);
		dSRT = std(tmpSRT);
		mv   = mean(tmpmv);
		dmv  = std(tmpmv);
	end
	
	fid = fopen(target_file,'w');
	fprintf(fid,'%5s\t %5s\t %5s\t %5s\n','SRT','dSRT','mv','dmv');
	fprintf(fid,'%2.2f\t %2.2f\t %1.5f\t %1.5f\n',SRT,dSRT,mv,dmv);
	fclose(fid);

end

% for a looped index
function out = loopedindex(in,maxind)
	out = mod(in-1,maxind)+1;
end

function [thewave,fs] = loadsentences(wavepath,idfiles)
	files = dir([wavepath filesep '*.wav']);
	out = [];
	for jj = 1:length(idfiles)
		[thewave,fs] = audioread([wavepath filesep files(idfiles(jj)).name]);
		out = [out; thewave];
	end
end

function [s,n] = res(s,n)
	n = n(loopedindex((1:length(s))+ceil(rand(1)*length(n)),length(n)));
end

function signal = crossfade_extend(signal, crossfade, len)
while (size(signal,1) < len)
  fade = repmat(linspace(0,1,crossfade).',1,size(signal,2));
  signal = [ ...
    signal(1:end-crossfade,:); ...
    signal(1:crossfade,:).* fade + ...
    signal(end-crossfade+1:end,:) .* (1-fade); ...
    signal(crossfade:end,:) ...
    ];
end
if size(signal,1) > len
	signal = signal(1:ceil(len),:);
end
end
