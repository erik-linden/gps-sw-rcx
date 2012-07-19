function ION = GEN_ION_LUT(file)

ION = {};

fid=fopen(file); %open file
if fid < 1
    warning('RAIN:NavData','Could not open file ''%s''.', sFileName);
    return;
end

% Find where the DCBs begins.
tline = fgetl(fid);
text = sscanf(tline,'%c',24);
while ~strcmp(text,'DIFFERENTIAL CODE BIASES')
    tline = fgetl(fid);
    text = sscanf(tline,'%c',24);
end

% Store the DCBs.
dcb = zeros(32,1);
for n = 1:32
    tline = fgetl(fid);
    if regexp(tline(4:6),'G[0-9][0-9]')
        ind = str2double(tline(5:6));
        tline(1:6) = [];
        dcb(ind) = sscanf(tline, '%f',1)*1e-9;
    end
end

% Find where the map begins.
tline = fgetl(fid);
tline(1:60) = [];
text = sscanf(tline,'%c',16);
while ~strcmp(text,'START OF TEC MAP')
    tline = fgetl(fid);
    tline(1:60) = [];
    text = sscanf(tline,'%c',16);
end

% Store the map.
time = zeros(13,1);
tec = zeros(73,71,13);
for tN = 1:13
    tline = fgetl(fid);
    ymdhms = sscanf(tline,'%f',6);
    gps_time = ymdhms_to_gps(ymdhms);
    time(tN) = gps_time(2);
    
    for latN = 1:71
        tline = fgetl(fid);
        for n = 0:3
            tline = fgetl(fid);
            values = sscanf(tline, '%f',16);
            tec(1+n*16:16+n*16,latN,tN) = values;
        end
        tline = fgetl(fid);
        values = sscanf(tline, '%f',9);
        tec(1+4*16:9+4*16,latN,tN) = values;
    end
    tline = fgetl(fid);
    tline = fgetl(fid);
end
fclose(fid);

% Create a structure.
ION.DCB = dcb;
ION.TEC = tec*1e-1; % EXPONENT set to -1
ION.Time = time';
ION.Height = 450e3; % HEIGHT set to 450 km
ION.LatGrid = 87.5:-2.5:-87.5;
ION.LonGrid = -180:5:180;