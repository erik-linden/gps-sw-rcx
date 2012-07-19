function [ymdhms] = gps_to_ymdhms(gps_time)
% [ymdhms] = gps_to_ymdhms(gps_time)
%
% Convert ymdhms time format to gps-week/seconds-of-week time format.
% 
% Input:
%   gps_time:     [gps week, second of week]
%
%   gps_week    - GPS week
%   sec_of_week - Second of the week
% 
% Output:
%   ymdhms:       [year, month, mday, hour, minute, second]
% 
%   year        - Year
%   month       - Number of the month
%   mday        - Day of the month
%   hour        - Hour of the day
%   minutes     - Minutes of the hour
%   seconds     - Seconds of the hour

JAN61980 = 44244;
JAN11901 = 15385;
SEC_PER_DAY = 86400;

gps_week = gps_time(1);
sec_of_week = gps_time(2);

month_day = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365;...
             0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366];

mjd = gps_week*7 + floor(sec_of_week/SEC_PER_DAY) + JAN61980;
fmjd = mod(sec_of_week, SEC_PER_DAY)/SEC_PER_DAY;

days_fr_jan1_1901 = mjd - JAN11901;
num_four_yrs = floor(days_fr_jan1_1901/1461);
years_so_far = 1901 + 4*num_four_yrs;
days_left = days_fr_jan1_1901 - 1461*num_four_yrs;
delta_yrs = floor(days_left/365) - floor(days_left/1460);

year = years_so_far + delta_yrs;
yday = days_left - 365*delta_yrs + 1;
hour = floor(fmjd*24);
minute = floor(fmjd*1440 - hour*60);
second = fmjd*86400 - hour*3600 - minute*60;

leap = (mod(year,4) == 0);
guess = floor(yday*0.032);
more = ((yday - month_day(1+leap,guess+2))  > 0);
month = guess + more + 1;
mday = yday - month_day(1+leap,guess+more+1);

ymdhms = [year, month, mday, hour, minute, second];