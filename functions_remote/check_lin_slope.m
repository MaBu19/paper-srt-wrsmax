function [r2_max,r2_idx] = check_lin_slope(levels,sii_tmp)

for n = 2:length(sii_tmp)-1
    [r, ~, ~, ~ ] = corrcoef(levels(n-1:n+1),sii_tmp(n-1:n+1));
    r2(n) = r(1,2).^2;
end
[r2_max,r2_idx] = max(r2);