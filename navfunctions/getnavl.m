function line = getnavl(fid)
% line = getnavl(fid) function to read nav file line and fix it up so that
% it can be read in Matlab (replace D with E and slit fields)
line = fgetl(fid);
if ischar(line)
    line = strrep(line,'D','E');
    if length(line) >= 79
        line = [line(1:22) ' ' line(23:41) ' ' line(42:60) ' ' line(61:79)]; 
    elseif length(line) >= 60
        line = [line(1:22) ' ' line(23:41) ' ' line(42:60) ' ' line(42:60)];
    elseif length(line) >= 41
        line = [line(1:22) ' ' line(23:41)];
    end
end