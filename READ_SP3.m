function S = READ_SP3(files)
% Reads a number of SP3 files and returns the results as a struc. If a
% satellite is not present, it will have all zeros in the data fields.
% Files must be in order.
%
% Input
%         files:  {'file1', 'file2'...}
%
% Output
% S =
%       GPSTime:  [1xEpocs]     Total seconds elapsed in the GPS time system
%             X:  [32xEpocs]    SV coordinates
%             Y:  [32xEpocs]
%             Z:  [32xEpocs]
%          dClk:  [32xEpocs]    SV clock offset

S = {};
weekOffset = 0;
i = 1; %counter of epochs for which the satellite coordinates are listed in PE
for npe = 1:length(files)
    %read header
    first ='1'; % variable for testing end of  header
    sFileName = files{npe};
    fid=fopen(sFileName); %open file
    if fid < 1
        warning('RAIN:NavData','Could not open file ''%s''.', sFileName);
        return;
    end
    tline = fgetl(fid);%skip 2 lines
    tline = fgetl(fid);
    tline = fgetl(fid);
    tline(1) = [];  %delete + character
    NSat = sscanf(tline, '%i',1); %read number of satellites
    while first ~= '*'
        tline = fgetl(fid);
        first = tline(1);
    end
    first ='1'; % variable for testing the line with time
    while feof(fid) == 0 && tline(1) ~= 'E';
        tline(1) = [];  %delete * character
        time = sscanf(tline,'%f',6);
        gps_time = ymdhms_to_gps(time);
        week = gps_time(1);
        seconds = gps_time(2);
        
        %Make sure data is in order
        if i==1
            firstWeek = week;
            S.GPSWeek = firstWeek;
        else
            dT = seconds - S.GPSTime(i-1);
            if (dT==(-604800+900))&&(week>(firstWeek+weekOffset))
                weekOffset = weekOffset + 1;
            elseif dT ~= 900
                warning(['Time increment between to two consecutive ephocs'...
                    ' was %d minutes, must be 15 minutes.'...
                    '\nMake sure input files are in order.'],dT/60)
            end
        end
        
        S.GPSTime(i) = 604800*weekOffset + seconds;
        
        tline = fgetl(fid);
        while first ~= '*'
            
            if tline(1) ~= 'P' % The code cant handle complicated SP3 files
                error('Unknown symbol "%s" encounterd in file %s, check format', tline(1), sFileName)
            end
            
            tline(1:2) = []; %delete PG
            prn = sscanf(tline(1:2), '%d'); %satellite number
            data = sscanf(tline,'%f',5);
            
            S.X(prn,i) = data(2)*1000;   %convert to [m]
            S.Y(prn,i) = data(3)*1000;
            S.Z(prn,i) = data(4)*1000;
            S.dClk(prn,i) = data(5)/1e6; %convert to [s]
            tline = fgetl(fid);
            if tline(1) == 'E'
                break;
            end
            first = tline(1);
        end
        i = i+1;
        first ='1';
    end
    fclose(fid);
end