function [levels,level_add] = add_levels(levels,sii_vec,step)

if step == 1
    dlev = [-2 2 4 6];
elseif step > 1
    dlev = [-1 1 2 3];
end


sii_diff = diff(sii_vec(levels));
[~, sii_diff_max] = max(sii_diff);   

if all(sii_diff==0) % all SII values the same -> not feasible to calculate slope
    level_add = [];
else
    level_add = levels(sii_diff_max)+dlev;
    level_add(level_add>27) = [];
    level_add(level_add<1) = [];
end

levels = sort(unique([levels, level_add])); % all current levels


