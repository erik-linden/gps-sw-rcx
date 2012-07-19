function [chip_rate_hist cst_err_hist cst_hist lock_hist magnitude_hist phi_err_hist...
          phi_if_hist rmse_hist tracking_hist w_df_err_hist w_df_hist] = ...
          LOAD_TRACK(sv)
% [chip_rate_hist cst_err_hist cst_hist lock_hist magnitude_hist phi_err_hist...
%   phi_if_hist rmse_hist tracking_hist w_df_err_hist w_df_hist] = ...
%   LOAD_TRACK(sv)

global TRACK_DIRECTORY COH_INT_TIME

S = load(sprintf('%stracking_hist_%d',TRACK_DIRECTORY,sv));

chip_rate_hist  = S.chip_rate_hist;
cst_err_hist    = S.cst_err_hist;
cst_hist        = S.cst_hist;
lock_hist       = S.lock_hist;
magnitude_hist  = S.magnitude_hist;
phi_err_hist    = S.phi_err_hist;
phi_if_hist     = S.phi_if_hist;
rmse_hist       = S.rmse_hist;
tracking_hist   = S.tracking_hist;
w_df_err_hist   = S.w_df_err_hist;
w_df_hist       = S.w_df_hist;

COH_INT_TIME = S.COH_INT_TIME;