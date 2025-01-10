function [percmat,xbins] = calc_percentiles_sections(xdata,ydata,xbins,filter_crit,percvec)
stepsize = (xbins(2)-xbins(1))/2; % assumes equidistant bins

for s = 1:length(xbins)
    percmat(s,:) = prctile(ydata(filter_crit & xdata > xbins(s)-stepsize & xdata <= xbins(s)+stepsize),percvec);
end 


