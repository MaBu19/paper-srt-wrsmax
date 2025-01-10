% function amplifies input signal by an amount specified in dB
%
% USAGE         out = ampSig(data, gain_dB)
%
% INPUT         data            signal to be amplified
%               gain_dB         scalar with gain in decibels (can be negative for attenuation)
%
% OUTPUT        out             amplified / attenuated signal
%
% Jan Rennies, 2006-03, rns@idmt.fraunhofer.de
function out = ampSig(data, gain_dB)

	out = data *10^(gain_dB/20);

end
