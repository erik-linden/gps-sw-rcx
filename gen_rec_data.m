% Script to generate dummy receiver data files. Use as reference
%  Output
%     Rec =
%          Time: [1xN] GPS system time at the start of the recording 
%       GPSTime: [1xN] GPS system time at data points
%           Pos: [3xN] Position at data points
%           Vel: [3xN] Velocity at data points
%          dFre: [1xN] Frequency offset at data points
%         Alpha: [1x1] System noise parameter
%           Eph: [Struct] Eph. structure
%       EphType: [1x1] Eph. type (0-3)
%           Ion: [Struct] Ion structure
%       IonType: [1x1] Ion type (0-2)


% START OF USER INPUTS
% -------------------------------------------------------------------

% Time at start of measurement, in year, month, day of mont, hour, minute,
% second. Use Swedish time.
year    = 2012;
month   = 04;
mday    = 27;
hour    = 15;
minute  = 38;
second  = 00;

% The location of the measurement, lon, lat and alt.
pos = lla2ecef([59.34967, 18.06910, 20])';% <- SPP
% pos = lla2ecef([59.3930089, 15.9251639, 20])';% <- ARBOGA
% pos = lla2ecef([48, 11, 20])';% <- DLR

% LO offset, 500 Hz has been common.
loOffset = 600;

% Which file to use when computing noise parameter.
nFile = 15;

% END OF USER INPUTS
% -------------------------------------------------------------------

% Convert to UTC. This will usually compensate for local SWE time.
if (month>3 && month<10) || (month==3 && mday>=25) || (month==10 && mday<=30)
    hour = hour-2;
else
    hour = hour-1;
end

% GPS time is UTC + 15 leap seconds.
second = second + 15;

% Compund date into a vector.
datevec = [year, month, mday, hour, minute, second];

% Compute the start time.
gpsTime = ymdhms_to_gps(datevec);
week = gpsTime(1);
startTime = gpsTime(2);

% Derive date information.
day = floor(startTime/(3600*24));
dayOfYear = floor(datenum(datevec)-datenum([year 0 0 0 0 0]));

% Print date information.
fprintf('\n..System time..\n')
fprintf('Week: %i, Day: %i\n',week,day)
fprintf('Day of year: %i\n',dayOfYear)

% Velocity defaults to zero.
vel = [0,0,0]';

% Where to put the endpoint.
length_of_data = 10*60;

% Compund the data into a structured array.
Rec.Time = [0, length_of_data];
Rec.GPSTime = [0, length_of_data]+startTime;
Rec.GPSWeek = week;
Rec.Pos = [pos, pos];
Rec.Vel = [vel, vel];
Rec.dFre = [loOffset, loOffset];
Rec.Alpha = 5e3;    % Typical value

%Download available data files.
fprintf('\n..Downloading data..\n')
[C, ephType, I, ionType] = DOWNLOAD_IGS(week,day,year,dayOfYear);

% Store the data.
Rec.Eph = C;
Rec.EphType = ephType;
Rec.Ion = I;
Rec.IonType = ionType;

% Set globals for ESTIMATE_DOPPLER
clear global;
global INITIALIZE RECEIVER_FILE
INITIALIZE = 1;
RECEIVER_FILE = 'rec';

% Save the rec-file to make it avaliable to ESTIMATE_DOPPLER.
save(RECEIVER_FILE,'Rec')

fprintf('\n..Estimating noise..\n')

if ~isempty(C)
    % Find the elevation of all SVs.
    elevArr = zeros(1,32);
    for n = 1:32
        [df, elev] = ESTIMATE_DOPPLER(0, n);
        elevArr(n) = elev;
    end

    % Sort them by elevation.
    [sortedElev sortedSv] = sort(elevArr);

    % Pick out those below 20 deg, but 10 at most.
    indLast = min(find(sortedElev<-20,1,'last'),10);
    
    fprintf('Using the %d lowest SV to compute alpha.\n',indLast)
else
    % There is no ephemeris data, so pick 10 SVs at random.
    sortedSv = randperm(32);
    indLast = 10;
    
    fprintf('Using 10 random SV to compute alpha.\n')
end

% Call CONSTANTS to get the file name.
clear global;
global FILE_NAME %#ok<NUSED>
CONSTANTS

% Estimate alpha...
alpha = estimate_noise_alpha(sortedSv(1:indLast),1,nFile); %nPrn = 2,

% ...and store it.
Rec.Alpha = alpha;

% Save the receiver file.
save(RECEIVER_FILE,'Rec')

% Clean up.
clear all