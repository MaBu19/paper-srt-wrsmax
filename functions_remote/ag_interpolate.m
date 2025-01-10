function [T_out] = ag_interpolate(T_in,ag_freqs)
ag_in = T_in{:,ag_freqs};
ag_out = ag_in;

num_freq = size(ag_in,2);
idx_del = false(size(ag_in,1),1);
for p = 1:size(ag_in,1)
    idx = find(isnan(ag_in(p,:)));
    if length(idx) >= num_freq-1
        idx_del(p)=true; % collect indices to delete (too few points)
    end
    if ~isempty(idx)
        if length(idx) < num_freq-1
            ag_i = interp1(1:num_freq,ag_in(p,:),idx,'spline');
            ag_out(p,idx) = max(min(ag_i,120),-10); % to avoid values outside audiometer limit range
        end
    end
end

ag_out = round(2*ag_out/10,0)*5;
T_in(:,ag_freqs) = table(ag_out(:,1),ag_out(:,2),ag_out(:,3),ag_out(:,4),ag_out(:,5),ag_out(:,6),ag_out(:,7),ag_out(:,8),ag_out(:,9),'VariableNames',ag_freqs);
T_in(idx_del,:) = [];
T_out = T_in;