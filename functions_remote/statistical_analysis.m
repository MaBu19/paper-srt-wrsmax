function [r2,p,bias,rms,rl2,ru2] = statistical_analysis(data_predicted, data_exp)
nanmask = isnan(data_predicted) | isnan(data_exp);
data_predicted = data_predicted(~nanmask);
data_exp = data_exp(~nanmask);
[r p rl ru] = corrcoef(data_predicted,data_exp);
r2 = r(1,2).^2;
p = p(1,2);
rl2 = sign(rl(1,2))*rl(1,2).^2;
ru2 = sign(ru(1,2))*ru(1,2).^2;
bias = fminsearch(@(x) my_costfun(data_predicted,data_exp,x),1);
rms = sqrt(mean((data_predicted-data_exp).^2));
end

function out = my_costfun(xin,yin,x)
yfun = x + 1*xin;
out = sqrt(mean(yin-yfun).^2);
end
