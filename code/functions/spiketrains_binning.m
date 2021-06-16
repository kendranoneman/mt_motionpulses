function [N_tot,N_b20] = spiketrains_binning(good_vals,ind_range)
%% Binning spikes into 1ms timebins
maxTime = 2851;
N_tot = [];
for l = 1:size(good_vals,1)
    ca = good_vals(l); ca2 = [ca{:}]';
    b20 = ca2(ca2(:) >= ind_range(1) & ca2(:) < ind_range(end));
    edges = 0 : 1 : (maxTime + 1);
    edges2 = ind_range(1) : 10 : ind_range(end);
    [N,edges] = histcounts(ca2,edges);
    [N2,edges2] = histcounts(b20,edges2);
    N_tot(l,:) = N;
    N_b20(l,:) = N2;
end
end