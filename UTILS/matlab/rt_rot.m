function [radial,trans] = rt_rot(north,east,cha)

% For rotating the waveform to RTZ 
% BAZ must be present in the header
%
% USAGE: [radial,trans] = rt_rot(north,east,'B')
% Input: north and east variables are read into MATLAB by rsac (n x3)

% check lengths of the two components
% note: if different by one point, then cut the last point
nlen = length(north);
elen = length(east);
if nlen ~= elen
    if abs(nlen-elen)==1
        disp('lengths differ by one point only -- cut the last point');
        if nlen > elen
            % north longer than east
            north(end,:) = [];                          % remove the last point
            north=ch(north, 'NPTS', elen);              % elen = nlen-1
            Enorth = lh(north,'E') - lh(north,'DELTA'); % modified end time
            north=ch(north,'E',Enorth);
        else
            % east longer than north
            east(end,:) = [];                           % remove the last point
            east=ch(east, 'NPTS', nlen);                % nlen = elen-1
            Eeast = lh(east,'E') - lh(east,'DELTA');    % modified end time
            east=ch(east,'E',Eeast);
        end
    else
        nlen, elen
        warning('lengths do not match');
        radial=0;trans=0;
        return
    end
    nlen = length(north);
    elen = length(east);
end

% To check that both are of same event (startTime and duration can be included)
if (lh(north,'BAZ') ~= lh(east,'BAZ'))
    error('Components have different backazimuth: rotation not done');
end

if (lh(north,'CMPINC') ~= 90 || lh(east,'CMPINC') ~= 90)
    error('Components are not horizontal');
end

% No rotation will be performed if CMPAZ information is missing
if (lh(north,'CMPAZ') == -12345) || (lh(east,'CMPAZ') == -12345)
    %disp(sprintf('Component Azimuth missing: considering them N and E aligned'));
    %north = ch(north,'CMPAZ',0);
    %east = ch(east,'CMPAZ',90);
    error('Component Azimuth missing');
else
    phi = abs(lh(north,'CMPAZ') - lh(east,'CMPAZ'));
	dif = mod(phi,90);
	if phi < 90 || (phi>180 && phi<270)
		dif = (90-dif);
	end
	if (phi ~= 90 && phi ~=270)
		disp('N and E are not orthogonal');
		if ((phi >= 86 && phi <= 94) || (phi >=266 && phi <= 274))
			disp(sprintf('N and E deviate from orthogonality by %f degree (maximum allowed 4 degrees)',dif));
        else
            warning('non orthogonal pair');
            radial=0;trans=0;
        return
		end
	end

% Checking if it is left-handed system or right
% LPSPOL TRUE if station components have a positive polarity (left-hand rule)
	if phi == 270
		if lh(east,'CMPAZ')>lh(north,'CMPAZ')
			east(:,2) = east(:,2)*-1;
			disp(sprintf('right-hand system'));
		end
	end
end
%==========================================================================
baz = lh(north,'BAZ');
cmpaz_i = lh(north,'CMPAZ');	% initial component azimuth
cmpaz_l = wrapTo360(baz+180);	% CMPAZ of radial component after rotation
phi = cmpaz_l - cmpaz_i;

% Estimating rotation angle
if (cmpaz_l > cmpaz_i)
	rot_angle = phi;
elseif (cmpaz_l < cmpaz_i)
	rot_angle = 360+phi;
else
	radial = north;
	trans = east;
	return
end

%-----Rotation-----
cphi = cosd(rot_angle);
sphi = sind(rot_angle);

n = north(:,2);
e = east(:,2);

r =  cphi*n + sphi*e;
t = -sphi*n + cphi*e;
%------------------

radial(:,1) = north(:,1);
radial(:,2) = r;
radial(:,3) = north(:,3);

trans(:,1) = east(:,1);
trans(:,2) = t;
trans(:,3) = east(:,3);

% Modifying header information
radial = ch(radial,'CMPAZ',cmpaz_l,'KCMPNM',[cha 'HR']);
trans = ch(trans,'CMPAZ',wrapTo360(cmpaz_l+90+dif),'KCMPNM',[cha 'HT']);

%==========================================================================
