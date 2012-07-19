function delta_w_df = PLL_2ND(cno, phi_err_k, phi_err_km)
% delta_w_df = PLL_2ND(cno, phi_err_k, phi_err_km)

global T COH_INT_TIME D1FDT1

% see http://nfudee.nfu.edu.tw/ezfiles/43/1043/img/327/ic5.pdf
Bn = ((9*cno^2*D1FDT1^2*pi^2*COH_INT_TIME)/(8+16*cno*COH_INT_TIME))^(1/5);  %get optimal bandwitdh
wn = sqrt(2)*4*Bn/3;

C1 = wn*(T*wn/2+sqrt(2));
C2 = wn*(T*wn/2-sqrt(2));

delta_w_df = C1*phi_err_k + C2*phi_err_km;   %estimated change in doppler frequency observed