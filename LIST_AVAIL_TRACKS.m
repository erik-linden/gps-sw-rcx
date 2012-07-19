function sv = LIST_AVAIL_TRACKS()
% Returns a sorted list of SV numbers for available tracks in the 
% TRACK_DIRECTORY.

global TRACK_DIRECTORY

listing = dir(TRACK_DIRECTORY);

nFiles = length(listing);

sv = [];

for n = 1:nFiles
    if (~listing(n).isdir)
        [token]=regexp(listing(n).name,'tracking_hist_(\d+).mat', 'tokens');
        
        if ~isempty(token)
            svNum = str2num(char(token{1})); %#ok<ST2NM>
            
            if (svNum>0 && svNum<33)
                sv = [sv svNum]; %#ok<AGROW>
            end
        end
    end
end

sv = sort(sv);