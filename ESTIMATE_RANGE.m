function range = ESTIMATE_RANGE(t, prn)

nSv = length(prn);
range = zeros(1, nSv);

for n = 1:nSv
    % Get GPS system time from receiver data.
    gps_time = GET_REC_DATA(t, 'time');
    
    % Get receiver position.
    recXYZ = GET_REC_DATA(t, 'pos');
    
    % Get SV position.
    svXYZ = GET_EPH(gps_time, prn(n), 'pos');
    
    % Calculate range
    range(n) = sqrt(sum((svXYZ-recXYZ).^2));
end

