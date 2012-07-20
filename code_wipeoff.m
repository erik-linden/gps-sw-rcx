CONSTANTS_H
CONSTANTS

load('tracking_hist_1')
t0 = floor(nanmin(cst_hist));

in_sig = LOAD_GPS_DATA(RAW_FILE, 120);
out_sig = zeros(size(in_sig));

code = CACODEGN(1);
index = BIT_LOCK(phi_if_hist, 1);
bit = BIT_EXTRACT(phi_if_hist, index);

l=length(bit);
bitu = bit(floor(1+(0:l*20-1)/20));
bitu = [zeros(index-1,1); bitu];


for n = 1:950;
    cst = cst_hist(n)-t0+floor(COH_INT_TIME/2e-3)*1e-3;
    if ~isnan(cst)
        samp = round(cst*FS);
        cst_offset = samp/FS-cst;
        
        corr = DIGITIZE_CA(code,cst_offset,ONE_MSEC_SAM);
        
        ind = samp+(1:ONE_MSEC_SAM);
        out_sig(ind) = corr.*in_sig(ind)*(-1)^bitu(n+floor(COH_INT_TIME/2e-3));
%         disp(bitu(n+floor(COH_INT_TIME/2e-3)))
%         disp(phi_if_hist(n+floor(COH_INT_TIME/2e-3)))
    end
end