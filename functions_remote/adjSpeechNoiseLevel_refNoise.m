function [tmp_speech, tmp_noise] = adjSpeechNoiseLevel_refNoise(tmp_speech, tmp_noise, SNR,Ref_Noise_Level)
	%% computation of the level 
	% USAGE    [tmp_speech, tmp_noise] =
	% adjSpeechNoiseLevel_refNoise(tmp_speech, tmp_noise, SNR,Ref_Noise_Level);
	%
	% INPUT     tmp_speech:         speech signal
	%           tmp_noise:          noise signal
	%           SNR:                SNR 
	%           Ref_Noise_Level:    level for noise, calibrated to dB FS 
	%
	% OUTPUT    tmp_speech:     speech signal, calibrated to dB FS (corresponding to needed SNR)
	%           tmp_noise:    	noise signal, calibrated to dB FS (rms 65 dB)
	%
	% Eugen Albertin, 2010-2, albeen@idmt.fraunhofer.de
	%
	%

	% set noise to reference level
	noise_mean_level = mean(20*log10(rms(tmp_noise))); % sound(tmp_noise)
	for ii=1:size(tmp_noise,2)
		  tmp_noise(:,ii) = ampSig(tmp_noise(:,ii),Ref_Noise_Level-noise_mean_level);
	end
	clear ii;
	% some checks noise_level_checks = 20*log10(rms(tmp_noise));

	% set speech corresponding to SNR Level
	speech_mean_level = mean(20*log10(rms(tmp_speech))); % sound(tmp_speech), figure; plot(tmp_speech)
	for ii=1:size(tmp_speech,2)
		  tmp_speech(:,ii) = ampSig(tmp_speech(:,ii),Ref_Noise_Level-speech_mean_level+SNR);
	end
	clear ii;

	% some checks
	if 0
		  noise_level_checks = 20*log10(rms(tmp_noise));
		  speech_level_checks = 20*log10(rms(tmp_speech));
		  SNR2 = 10*log10(sum(tmp_speech(:,1).^2)./sum(tmp_noise(:,1).^2));
	end
end
