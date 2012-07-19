function [gps_time] = ymdhms_to_gps(ymdhms)
% [gps_time] = ymdhms_to_gps(ymdhms)
%
% Convert gps-week/seconds-of-week time format to ymdhms time format.
% 
% Input:
%   ymdhms:       [year, month, mday, hour, minute, second]
% 
%   year        - Year
%   month       - Number of the month
%   mday        - Day of the month
%   hour        - Hour of the day
%   minutes     - Minutes of the hour
%   seconds     - Seconds of the hour
%
% Output:
%   gps_time:     [gps week, second of week]
%
%   gps_week    - GPS week
%   sec_of_week - Second of the week

JAN61980 = 44244;
JAN11901 = 15385;
SEC_PER_DAY = 86400;

year = ymdhms(1);
month = ymdhms(2);
mday = ymdhms(3);
hour = ymdhms(4);
minute = ymdhms(5);
second = ymdhms(6);

month_day = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334;...
             0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];

leap = (mod(year,4) == 0);
yday = month_day(1+leap,month) + mday;

mjd = (floor((year - 1901)/4))*1461 + mod((year - 1901),4)*365 + yday - 1 + JAN11901;
fmjd = ((second/60 + minute)/60 + hour)/24;

gps_week = floor((mjd - JAN61980)/7);
sec_of_week = ( (mjd - JAN61980) - gps_week*7 + fmjd )*SEC_PER_DAY;

gps_time = [gps_week, sec_of_week];