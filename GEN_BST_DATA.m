function GEN_BST_DATA(prn)
% Creates a file with the times new bits starts for each SV. Useful when
% coherent integration time is more than one millisecond.

bst = NaN(32,1);
nSv = length(prn);

for n = 1:nSv
    load(sprintf('tracking_hist_%d',prn(n)));

    I = BIT_LOCK(phi_if_hist,prn(n));
    
    x = 1:length(cst_hist);
    indGood = ~isnan(cst_hist);
    
    bst(prn(n)) = interp1(x(indGood),cst_hist(indGood),I,'linear','extrap');
end

save('bst.mat','bst')