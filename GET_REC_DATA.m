function var = GET_REC_DATA(tRel, string)
% var = GET_REC_DATA(tRel, string)
%
% Returns interpolated receiver data.
% 
% Input:
%             tRel:    time reletive start of recording
%           string:    string specifying which data to return,
%                      'time' gives current GPS system time
%                      'pos'  gives position
%                      'vel'  gives velocity
%                      'dFre' gives LO frequency offset
% 
% Output:
%              var:     [1x3] for 'vel' and 'pos', scalar for 'time' and 'dFre'

global INITIALIZE RECEIVER_FILE
persistent Rec

% Make sure data file is loaded.
if INITIALIZE
    load(RECEIVER_FILE,'Rec');
end

switch string
    case 'time'
        var = quick_interp(Rec.Time, Rec.GPSTime, tRel,1);
    case 'pos'
        var = zeros(1,3);
        var(1) = quick_interp(Rec.Time, Rec.Pos(1,:), tRel);
        var(2) = quick_interp(Rec.Time, Rec.Pos(2,:), tRel);
        var(3) = quick_interp(Rec.Time, Rec.Pos(3,:), tRel);
    case 'vel'
        var = zeros(1,3);
        var(1) = quick_interp(Rec.Time, Rec.Vel(1,:), tRel);
        var(2) = quick_interp(Rec.Time, Rec.Vel(2,:), tRel);
        var(3) = quick_interp(Rec.Time, Rec.Vel(3,:), tRel);
    case 'dFre'
        var = quick_interp(Rec.Time, Rec.dFre(1,:), tRel);
    otherwise
        error('No receiver data for "%s"',string)
end