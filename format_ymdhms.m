function s = format_ymdhms(ymdhms, swe)
% s = format_ymdhms(ymdhms, Swe)
%
% Returns a string of formated ymdhms date.


if swe %Swedish local time
    ymdhms(4) = ymdhms(4) + 2;
    ymdhms(6) = ymdhms(6) - 15;
    type = 'UTC+2';
else %GPS system time
    type = 'GPS';
end

s = sprintf('%s Year:%d Month:%d Day:%d Hour:%d Minute:%d Second:%f',...
    type, ymdhms);