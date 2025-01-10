% maximum discrimination (in speech tests such as FMST)

function [SI_max,level_max,SI_min,level_min] = max_SI_discr(SI_vec,level_vec)

if size(SI_vec,2) == length(level_vec)
    % maximum
    [SI_max,level_idx] = max(SI_vec,[],2);

    level_max = level_vec(level_idx);
    level_max(isnan(SI_max)) = nan;

    % minimum
    [SI_min,level_idx_min] = min(SI_vec,[],2);

    level_min = level_vec(level_idx_min);
    level_min(isnan(SI_min)) = nan;

end
