function [w_df_expObs, w_df_expAct, delta_w_df_expObs, delta_w_df_expAct, dphidt] = GET_TRACKING_DATA(t, prn)
% [w_df_expObs, w_df_expAct, delta_w_df_expObs, delta_w_df_expAct, dphidt] = GET_TRACKING_DATA(t, prn)
% 
% Input:
%   t                 - time relative start of recording 
%   prn               - SV PRN
% 
% Output:
% 	w_df_expObs       - predicted apperent Doppler frequency
%   w_df_expAct       - predicted actual Doppler frequency
%   delta_w_df_expObs - predicted apperent change in Doppler frequency
%   delta_w_df_expAct - predicted actual change in Doppler frequency
%   dphidt            - predicted mean time derivative of Doppler phase

global T

w_df_1 = 2*pi*ESTIMATE_DOPPLER(t+T/2, prn);             %expected dopple frequency at middle of current code
w_df_2 = 2*pi*ESTIMATE_DOPPLER(t+T*3/2, prn);           %expected dopple frequency at middle of next code

w_dLO1 = 2*pi*GET_REC_DATA(t+T/2, 'dFre');              %expected receiver LO offset at middle of current code
w_dLO2 = 2*pi*GET_REC_DATA(t+T*3/2, 'dFre');            %expected receiver LO offset at middle of next code

w_df_expObs = w_df_1 - w_dLO1;                          %expected observed doppler frequency
w_df_expAct = w_df_1;                                   %expected actual doppler frequency

delta_w_df_expObs = (w_df_2-w_dLO2)-(w_df_1-w_dLO1);    %expected change in observed doppler
delta_w_df_expAct = (w_df_2-w_df_1);                    %expected change in actual doppler frequency
dphidt = delta_w_df_expAct/2;                           %expected mean phase change rate