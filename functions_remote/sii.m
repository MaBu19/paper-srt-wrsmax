function [SII]   = sii(speech_signal, noise_signal, fs ,f_audiogram, dBHL_audiogram)

	%% computation of the SII
	% USAGE     [SII results_SII]   = my_my_STI(tmp_speech, tmp_noise, fs,...
	%                                   f_audiogram, dBHL_audiogram)
	%
	% INPUT     speech_signal:     speech signal, calibrated to dB FS
	%           noise_signal:      noise signal, calibrated to dB FS
	%           fs:                sampling frequency
	%           f_audiogram:       frequencies of the measured audiogramm
	%           dBHL_audiogram:    audiogramm data to to corresponding f_audiogram
	%
	% OUTPUT    SII:      	SII
	%           results:    structure containing intermediate results
	%
	% Eugen Albertin, 2010-2, albeen@idmt.fraunhofer.de
	%

	% frequency filtering into critical bands, using design of octave-package
	bands = [ 100   200   150 ;
		        200   300   250 ;
		        300   400   350 ;
		        400   510   450 ;
		        510   630   570 ;
		        630   770   700 ;
		        770   920   840 ;
		        920  1080  1000 ;
		        1080  1270  1170;
		        1270  1480  1370;
		        1480  1720  1600;
		        1720  2000  1850;
		        2000  2320  2150;
		        2320  2700  2500;
		        2700  3150  2900;
		        3150  3700  3400;
		        3700  4400  4000;
		        4400  5300  4800;
		        5300  6400  5800;
		        6400  7700  7000;
		        7700  9500  8500];

	weight = 'spin';     % weighting for SPIN

	fc_bank = bands(:,3);
	flow  = bands(:,1); 
	fhigh = bands(:,2); 

	N = 3;                                          % filter order
	filt_sig = zeros(length(speech_signal),length(fc_bank));
	filt_noise = zeros(length(noise_signal),length(fc_bank));
	% filt_ir = zeros(size(in_ir));
	I_speech = zeros(size(fc_bank));
	I_noise = zeros(size(fc_bank));

	filt_sig   = Critical_band_filterbank_o3d(speech_signal,fs);
	filt_noise = Critical_band_filterbank_o3d(noise_signal,fs);
	I_speech = rms(filt_sig.^2,1)';
	I_noise = rms(filt_noise.^2,1)';
	
	I_speech_dB = 10*log10(I_speech);
	I_noise(I_noise==0) = realmin;
	I_noise_dB = 10*log10(I_noise);
	% if  strmatch('none', noise_type,'exact')
	%     I_noise_dB = -100*ones(size(I_noise_dB));
	% end

	% compute power spectral density (PSD) from filter bandwidths and
	% intensities, and FF-2-Eardrum transformation
	fSpeech = pre_SII(I_speech_dB);
	fNoise = pre_SII(I_noise_dB);



	% interpolate audiogram data to get same vector length as speech spectrum
	% level:
	fHL = interp1(f_audiogram, dBHL_audiogram,fc_bank,'linear','extrap');

	SII = my_SII(fSpeech,fNoise,fHL,weight,0);

end
