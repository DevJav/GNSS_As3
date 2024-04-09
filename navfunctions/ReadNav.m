function [nav, leapsec, cnt, ion] = ReadNav(navfile)
% [nav leapsec cnt] = ReadNav(navfile) reads a rinex broadcast ephemeris file
% and return results in structure nav, number of leap seconds in leapsec
% and the number of entries in cnt
%
% The returned nav structure has the following elements
% .prn     - prn is the PRN number of satellite
% .date    - yyyy mm dd hh mm sec is the epoch to which the clock parameters
%                       follow apply
% .afc     - af0, af1, af2 are the clock offset (sec), rate (sec/sec) and
%                       acceleration (sec/sec^2).
% .a_f0     - af0
% .a_f1
% .a_f2
% .adoe    - adoe       age of ephemeris entry (sec) i.e, how long ago was 
%                       it uploaded
% .C_rs
% .C_rc    - crs,crc    radius corrections sin and cos  (m)
% .deltan  - dn         correction to mean motion (radian/sec)
% .M_0     - M0         Mean anomaly (radians)
% .C_uc
% .C_us    - cuc,cus    correction to argument in latitude sin and cos(rad)
% .e       - ecc        eccentricity
% .sqrtA   - art        square root of semi-major axis (m^.5)
% .t_oe    - toe        time of ephemeris (same as date+sec)
% .C_ic
% .C_is    - cic,cis    corrections to inclination sin and cos (rad)
% .omega_0 - om0        longitude of the ascending node (rad) (Capital omega)
% .i_0     - i0         inclination (rad)
% .omega   - w          argument of perigee (rad) (lower case omega)
% .omegaDot - omd        time derivative of longitude of the ascending node (rad/sec)
% .iDot    - idt        time derivative of inclination (rad/sec)
% .cflg12
% .pflg12  - cflg12, pflg12 Flags (whose meaning is not clear)
% .weekno  - weekno     GPS Week Number
% .svacc   - svaccuracy range accuracy (m)
% .svheath - svhealth   satellite health flag
% .T_GD    - tgd        Group delay L2 bias (word 7 subframe 1)
% .aodc    - aodc       age of clock parameter upload (sec)
% .t_oc    - transmit   Transmission time seconds of GPS week

% Script from: http://www-gpsg.mit.edu/~tah/12.540/ReadNav.m
% Ionosphere parameter extraction added by Daniel Olesen, DTU Space


nav = []; leapsec = 0;

% Try to open the navfile
fid = fopen(navfile);
if fid < 1 
    fprintf('Unable to open %s, nav structure returned empty\n',navfile);
    return
end

% Make sure this is a nav file
line = fgetl(fid);
ind = findstr(line,'RINEX VERSION');
if ind > 0 
    rxver = sscanf(line,'%f',1);
    fprintf('Rinex version %f found in %s\n',rxver, navfile);
else
    fprintf('%s in valid rinex file\nFirst line is:\n%s\n',navfile,line);
    return
end

if nargout > 3 % If  Ionospheric parameters are requested
    ion = struct;
    % skip the next two lines of header
%    fgetl(fid); 
    fgetl(fid);
    line = fgetl(fid); % Read alpha parameters
    line = strrep(line,'D','E');
    val= sscanf(line,'%f');
    ion.alpha0 = val(1); ion.alpha1 = val(2);
    ion.alpha2 = val(3); ion.alpha3 = val(4);
    
    line = fgetl(fid); % Read beta parameters
    line = strrep(line,'D','E');
    val = sscanf(line,'%f');
    ion.beta0 = val(1); ion.beta1 = val(2);
    ion.beta2 = val(3); ion.beta3 = val(4);
    
    
end

% Now skip rest of header
stillhead = 1;
while stillhead
    line = fgetl(fid);
    ind = findstr(line,'LEAP SECONDS');
    if ind > 0
        leapsec = sscanf(line,'%d',1);
    end
    ind = findstr(line,'END OF HEADER');
    if ind > 0 , stillhead = 0; end
end

% Now we are ready to read the file
cnt = 0;  % cnt will count number of entries
OK = 1 ;  % Set OK true; when EOF reach set to false.
while OK
    line = getnavl(fid); % Reads and fixes line.
    if ~ischar(line) 
        OK = 0;  % Return was empty and so EOF reached
    else
        vals = sscanf(line,'%f');
        % Save values;
        % prn yy mm dd hh mm sec     af0          af1          af2
        prn = vals(1);
        date = vals(2:7); % Year mon day hr min sec
        if date(1) < 50 
            date(1) = 2000+date(1);
        else
            date(1) = 1900+date(1);
        end
        afc = vals(8:10);
        % OK we are done with prn and date line, read the next block of lines
        line = getnavl(fid);
        vals = sscanf(line,'%f');
        % aode               crs         dn           M0
        aode = vals(1); crs = vals(2); dn = vals(3); M0 = vals(4);
        % Next line
        line = getnavl(fid);
        vals = sscanf(line,'%f');
        % cuc                ecc         cus          art
        cuc = vals(1); ecc = vals(2); cus = vals(3); art = vals(4);

        line = getnavl(fid);
        vals = sscanf(line,'%f');
        % toe                cic         om0          cis
        toe = vals(1); cic = vals(2); om0 = vals(3); cis = vals(4);
        line = getnavl(fid);
        vals = sscanf(line,'%f');
        % i0                 crc         w            omd
        i0 = vals(1); crc = vals(2); w = vals(3); omd = vals(4);
        line = getnavl(fid);
        vals = sscanf(line,'%f');
        % idt                cflg12      weekno       pflg12
        idt = vals(1); cflg12 = vals(2); weekno = vals(3); pflg12 = vals(4);

        line = getnavl(fid);
        vals = sscanf(line,'%f');
        % svaccuracy         svhealth    tgd          aodc
        svacc = vals(1); svhealth = vals(2); tgd = vals(3); aodc = vals(4);
        line = getnavl(fid);
        vals = sscanf(line,'%f');
        % transmit           spare       spare        spare
        transmit = vals(1);
        % Now construct the structure
        
        % Modified by Daniel
        
        a_f2 = afc(3);
        a_f1 = afc(2);
        a_f0 = afc(1);
        %
        
        str = struct('prn',prn,'date',date,'afc',afc, ...
            'a_f2',a_f2,'a_f1',a_f1,'a_f0',a_f0, ...
            'aode',aode,'C_rs',crs,'deltan',dn,'M_0',M0, ...
            'C_uc',cuc,'e',ecc,'C_us',cus,'sqrtA',art, ...
            't_oe',toe,'C_ic',cic','omega_0',om0,'C_is',cis, ...
            'i_0',i0,'C_rc',crc,'omega',w,'omegaDot',omd, ...
            'iDot',idt,'cflg12',cflg12,'weekno',weekno,'pflg12',pflg12, ...
            'svacc',svacc,'svhealth',svhealth,'T_GD',tgd,'aodc',aodc, ...
            't_oc',transmit);
        % Now add to nav array
        cnt = cnt + 1; nav = [nav str];
    end
end
        

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



        
        
       
    

        
     