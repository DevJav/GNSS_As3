function [eph,ion] = importGPSephemeris(navFile, tow)
%
% Import ephemeris-data from Nav file for epoch closest to tow (within one
% hour)
%
% 2016, Daniel Olesen DTU Space

% import broadcast ephemerides from Rinex NAV file

if nargout > 1
[nav,leapsec,cnt,ion] = ReadNav(navFile);
else
[nav,leapsec,cnt] = ReadNav(navFile);
end    
% read out the valid set of ephmerides based on mission time
% method should be altered for very long missions

nav = rmfield(nav,'date');
eph = struct;
%ion = struct; 

t_mission_start = tow;
for prn=1:32
    % find data for sv in struct
    idx = find([nav.prn]==prn);
    idx2 = find(abs([nav(idx).t_oe]-t_mission_start) < 3600);
    
    if (length(idx2)==1)
        s = sprintf('PRN%02d',prn);
        eph.(s) = nav(idx(idx2));
    end

    
end