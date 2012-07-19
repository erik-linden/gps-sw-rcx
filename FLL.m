function [delta_w_df, w_df_dot_kp] = FLL(w_df_dot_k, phi_act_k, phi_act_km)
% [delta_w_df, w_df_dot_kp] = FLL(w_df_dot_k, phi_act_k, phi_act_km)
%
% Frequency locked loop for tracking.
% 
% Input:
%     w_df_dot_k  - current frequency change rate [rad/s]
%     phi_act_k   - current phase angle relative local carrier
%     phi_act_kp  - previous phase angle relative local carrier
%     
% Output
%     delta_w_df  - frequency offset for local carrier
%     w_df_dot_kp - next frequency change rate [rad/s]

global T A_FLL B_FLL

% sin(theta) where theta is the phase change from km to k
rotation_angle = sin(diff(unwrap([phi_act_km phi_act_k])));

% Next frequency change rate is the current rate plus a constant times the
% rotation angle.
w_df_dot_kp = w_df_dot_k + A_FLL*rotation_angle;        

delta_w_df = w_df_dot_k*T + B_FLL*rotation_angle;