function GEN_BST_DATA(prn)
% Creates a file with the times new bits starts for each SV. Useful when
% coherent integration time is more than one millisecond.

global BIT_START_TIME_FILE

bst = NaN(32,1);
nSv = length(prn);

for n = 1:nSv
    [~, ~, cst , ~, ~, ~, phi_if , ~, ~, ~, ~] = LOAD_TRACK(prn(n));

    I = BIT_LOCK(phi_if,prn(n));
    
    x = 1:length(cst);
    indGood = ~isnan(cst);
    
    bst(prn(n)) = interp1(x(indGood),cst(indGood),I,'linear','extrap');
end

save(BIT_START_TIME_FILE,'bst')