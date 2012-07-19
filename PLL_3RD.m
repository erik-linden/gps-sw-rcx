function delta_w_df = PLL_3RD(snr, w_df_k, w_df_km, phi_err_k, phi_err_km, phi_err_km2)
% delta_w_df = PLL_3RD(...
% snr, w_df_k, w_df_km, phi_err_k, phi_err_km, phi_err_km2)

global T COH_INT_TIME D2FDT2

% see http://nfudee.nfu.edu.tw/ezfiles/43/1043/img/327/ic5.pdf
Bn = ((15625*snr^2*D2FDT2^2*pi^2*COH_INT_TIME)/(1458*(1+2*snr*COH_INT_TIME)))^(1/7);
wn = 1.2*Bn;

C1 = wn*(2+T*wn+T^2*wn^2/4);
C2 = wn*(-4+T^2*wn^2/2);
C3 = wn*(2-T*wn+T^2*wn^2/4);

delta_w_df = w_df_k - w_df_km + C1*phi_err_k + C2*phi_err_km + C3*phi_err_km2;