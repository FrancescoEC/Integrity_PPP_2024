function [trop_dry, trop_wet] = tropunb1(lla,t_gps,elev);

% [trop_dry, trop_wet] = tropunb1(lla,t_gps,elev);
%
% Computes the wet and dry troposphere delays using
% the University of New Brunswick (UNB1), altitude dependent, 
% troposphere model.
%
% The UNB1 model is a composite model that uses the explicit 
% forms of Saastamoinen's delay algorithms combined with 
% Niell's mapping functions. The surface and lapse rate parameters 
% are from the U.S. 1976 Standard Atmosphere. There will be a
% several centimetre bias at the poles and the equator unless 
% surface met data is used.  This model was developed at the
% University of New Brunswick under contract to Nav Canada.
%
% Input:
%   lla        - matrix of geodetic position (1x3 or nx3) [lat, lon, height]
%                 lat and lon are in rad, height is in m
%                 valid latitudes are -pi/2 -> pi/2
%                 valid longitudes are -pi -> 2*pi
%   t_gps      - GPS time vector for az/el pairs [GPS_week, GPS_sec] (nx2)
%                or [GPS_week GPS_sec rollover_flag] (nx3) for times prior
%                to Aug. 22, 1999.
%                 valid GPS_week values are 1-1024
%                 valid GPS_sec values are 0-604799
%   elev       - elevation angle to GPS satellites (rad) (nx1)
%                 Valid elevations are between -pi/2 and pi/2.
%						Elevations below .1 degree will have the a 
%        			delay mapped to .1 degree.  No mapping below
% 						zero degree elevation is modeled.
% Output:
%   trop_dry   - Dry troposphere component (m) (nx1)
%   trop_wet   - Wet troposphere component (m) (nx1)
%
% See also PSEUDO_R 

% Written in MATLAB by: Jimmy LaMance 4/12/99
% Copyright (c) 1999 by Constell, Inc.
%
% Original FORTRAN by: Paul Collins, University of New Brunswick
% UNB intellectual property rights reserved.
% For more information on UNB visit:
% http://www.unb.ca/GGE/Personnel/Langley/Langley.html

% Reference for UNB1 Troposphere Model: 
% Collins, J.P. and R.B. Langley. "A Tropospheric Delay Model for the User of
% the Wide Area Augmentation System." Final contract report prepared for Nav
% Canada, Department of Geodesy and Geomatics Engineering Technical Report
% No. 187, University of New Brunswick, Fredericton, N.B., Canada.
%
% Available in PDF format from
% http://gauss.gge.unb.ca/papers.pdf/waas.tropo.oct96.pdf

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
trop_dry=[]; 
trop_wet=[]; 

% Check the number of input arguments and issues a message if invalid
msg = nargchk(3,3,nargin);
if ~isempty(msg)
  fprintf('\n');
  dbstack
  fprintf('\n');
  fprintf('%s  See help on TROPUNB1 for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

estruct.func_name = 'TROPUNB1';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 't_gps';
estruct.variable(1).req_dim = [901 2; 901 3];
estruct.variable(1).var = t_gps;
estruct.variable(1).type = 'GPS_TIME';
  
estruct.variable(2).name = 'lla';
estruct.variable(2).req_dim = [901 3; 901 2; 1 3; 1 2];
estruct.variable(2).var = lla;
estruct.variable(2).type = 'LLA_RAD';
  
estruct.variable(3).name = 'elev';
estruct.variable(3).req_dim = [901 1];
estruct.variable(3).var = elev;
estruct.variable(3).type = 'ELEVATION_RAD';
  
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

% Set constants for T, P, and e. 
tempkelv = 288.15 - 0.00650 * lla(:,3);
I_zero = find(tempkelv < .1);
if ~isempty(I_zero)
   tempkelv(I_zero) = ones(size(tempkelv(I_zero)))*.1;
end

presmbar = 1013.25 *(tempkelv ./ 288.15).^5.256;
wvpmbar  = 11.691 *(tempkelv ./ 288.15).^21.024;
      
% Convert the elevation angle to degrees for future use
I_zero = find(elev < .0017);
if ~isempty(I_zero)
   elev(I_zero) = ones(size(elev(I_zero)))*.0017;
end
elev = elev * 180/pi;

% Compute the vertical delays for the wet and dry tropo
%   P,T,e equivalent to 324.8 N units at surface/msl
%   Uses Essen & Froome constants for consistency with Saastamoinen
[hzd,wad] = tzdds87(presmbar,tempkelv,wvpmbar,0.0065,3.0,lla);
     
% Compute the day of year from the current GPS time 
[ut, ls, day_of_year] = gps2utc(t_gps);   						% integer day
day_of_year = day_of_year + rem(t_gps(:,2),86400) ./ 86400;	% int + fractional day

% Compute the mapping function for the dry tropo
latdeg = lla(:,1)*180/pi;
hmf = nhmf2(day_of_year,latdeg,lla(:,3),elev);

% Compute the mapping function for the wet tropo
wmf = nwmf2(latdeg,elev);
     
% Compute the total troposphere delay
trop_dry = hzd.*hmf(:,1);
trop_wet = wad.*wmf(:,1);

% end of tropunb1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Built-in subroutines for UNB tropo model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tzd1, tzd2] = tzdds87(P, T, E, t_lapse, lambda, lla);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Comments from the original FORTRAN code.
%
%     FUNCTION: Computes the hydrostatic zenith delay using the Davis,
%               et.al. formulation of Saastamoinen's model and the wet
%               delay based on an updated 'Smith & Weintraub' expression
%               and incorporating the two-parameter formula first derived
%               by Saastamoinen and then Askne & Nordius.
%
%     REFERENCES: (1) Davis,J.L., T.A. Herring, I.I. Shapiro, A.E.E. Rogers
%                 and G. Elgered. (1985). "Geodesy by radio interferometry:
%                 Effects of atmospheric modelling errors on estimates of
%                 baseline length." Radio Science, Vol. 20, No. 6, pp. 1593
%                  - 1607.
%
%                 (2) Askne,J. and H. Nordius (1987). "Estimation of
%                 tropospheric delay for microwaves from surface weather
%                 data." Radio Science, Vol. 22, No. 3, pp. 379 - 386.
%
%                 (3) Saastamoinen, J. (1973). "Contributions to the theory
%                 of atmospheric refraction." In three parts. Bulletin
%                 Geodesique, No. 105, pp.279-298; No. 106, pp.383-397;
%                 No. 107, pp.13-34.
%
%     HISTORY: 
%               ORIGINAL CODE:    Paul Collins        19.05.96
%   
%     PARAMETERS AND VARIABLES:
%
%         INPUT:
%                P        - Total Atmospheric Pressure          [mbar]
%                T        - surface temperature                 [K]
%                E        - surface water vapor pressure        [mbar]
%                t_lapse   - temperature lapse rate              [K/m]
%                lambda   - the "lambda" parameter         
%                height       - height of the site                  [m]
%                latitude      - geodetic latitude                   [rad]
%
%         LOCAL:
%                rd       - molar gas constant of dry air       [J/(K.kg)]
%                k1,K3'   - Refractivity constants:
%                  k1     - constant for dry/moist gas          [K/mbar]
%                  K3'    - derived constant for water vapour   [K^2/mbar]
%
%                Values are from Thayer (1974) quoted in Davis.
%                k1                = 77.604 +/- 0.014
%                k3' = k3 + k2'.Tm = 382000 +/- 4000
%
%                c1       - constant for Rd, g0, k1
%                c2       - constant for Rd, k3'
%
%        OUTPUT:
%                tzd1     - hydrostatic zenith delay             [m]
%                tzd2     - wet zenith delay                     [m]
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get height and latitude from lla
height = lla(:,3);
latitude = lla(:,1);

% Set up constants
k1 = 77.604;
k3prime = 3.82d5;
g0 = 9.784;
rd = 287.054;
c1 = 1e-6 * rd * k1 / g0;
c2 = 1e-6 * rd * k3prime;

%-----------------------------------------------------------------------
%        compute the acceleration at the mass center
%        of a vertical column of the atmosphere
%-----------------------------------------------------------------------
dgref = 1.0 - 2.66d-03*cos(2.0 *latitude) - 2.8e-07 * height;

%-----------------------------------------------------------------------
%               compute zenith hydrostatic delay
%-----------------------------------------------------------------------
tzd1 = c1 ./ dgref .* P;

%-----------------------------------------------------------------------
%               compute zenith wet delay
%-----------------------------------------------------------------------
tzd2 = c2 ./ (g0 .* dgref *(lambda + 1.0) - rd * t_lapse) .* E ./ T;

% end of tzdds87

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin nwmf2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wmf = nwmf2(latitude, elev)
 
% new aen 930517 Routine to compute the new wmf2.0 mapping function which
%               depends only on latitude. 
 
% Define variables used 
%  a,b,c       - the a,b,and c coefficients in the continued fraction
%                form of Marini
%  beta        - intermediate term in calculation
%  gamma       - intermediate term in calculation
%  sine        - Sine of elevation angle
%  cose        - Cos of elevation angle
%  wmf(1)      - wet delay mapping function
%  wmf(2)      - d_wet_mapping_function/d_elevation
%  topcon      - Constant of top of mapping fuinction to ensure
%                that value is 1.0000 at zenith 
%  latitude   - latitude (degrees)
%  l          - absolute latitude
%  dl         - incremental latitude from last lat_wmf
%  elev       - elevation (degrees)
%  dl,da,db,dc  - used for interpolation

%  define parameters used for calculating coefficients.

lat_wmf = [15 30 45 60 75];

%  coefficients are from fits to raytraces of the standard atmospheres
%  for July for latitudes 15, 45, 60, and 75 degrees latitude and for 
%  January for 30 degrees latitude (930517).

abc_w2p0 = [5.8021897e-4,5.6794847e-4,5.8118019e-4,5.9727542e-4,6.1641693e-4; ...
     1.4275268e-3,1.5138625e-3,1.4572752e-3,1.5007428e-3,1.7599082e-3; ...
     4.3472961e-2,4.6729510e-2,4.3908931e-2,4.4626982e-2,5.4736038e-2]';

deg2rad = 3.14159265/180;

l = abs(latitude);

%  Coefficients for the continued fraction expansion for each latitude.
%  for latitudes less than 15 degrees:
I_15 = find(l <= lat_wmf(1));
if ~isempty(I_15)
   a(I_15) = abc_w2p0(1,1);
	b(I_15) = abc_w2p0(1,2);
	c(I_15) = abc_w2p0(1,3);
end % if ~isempty(I_15)

%  for latitudes between 15 and 75  degrees:
for i = 1:4
   I_mid = find(l > lat_wmf(i) & l <= lat_wmf(i+1));
   if ~isempty(I_mid)
      dl = (l(I_mid)-lat_wmf(i))/(lat_wmf(i+1)-lat_wmf(i));
		da = abc_w2p0(i+1,1)-abc_w2p0(i,1);
		a(I_mid) = abc_w2p0(i,1) + dl*da;
		db = abc_w2p0(i+1,2)-abc_w2p0(i,2);
		b(I_mid) =   abc_w2p0(i,2) + dl*db; 
		dc =   abc_w2p0(i+1,3)-abc_w2p0(i,3);
		c(I_mid) =   abc_w2p0(i,3) + dl*dc; 
	end % if ~isempty(I_mid)
end % for i = 1:4

%  for latitudes greater than 75 degrees:
I_high = find(l > lat_wmf(5));
if ~isempty(I_high)
   a(I_high) = abc_w2p0(5,1);
	b(I_high) = abc_w2p0(5,2);
	c(I_high) = abc_w2p0(5,3);
end % if

%  Now the coefficients exist; calculate the mapping function, wmf(1),
%      and the change of mapping function with elevation,
%      dwmf/d_el =wmf(2).
%  To calculate the delay-rate correction, d_tau/dt:
%      d_tau/dt = d_tau_zen/dt * wmf(1) + tau_zen * dwmf/d_el * d_el/dt 
sine = sin( elev * deg2rad);
cose = cos( elev * deg2rad);
beta = b' ./ (sine + c');
gamma = a' ./ (sine + beta);
topcon = (1.0 + a' ./ (1.0 + b' ./ (1.0 + c')));

wmf(:,1) = topcon ./ ( sine + gamma );

wmf(:,2) = -topcon ./ (sine + gamma).^2 .* ...
	( cose - a' ./ (sine + beta).^2 .* cose .* ...
	( 1.0 - b' ./ (sine + c').^2 ));

% end of nwmf2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin nhmf2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hmf = nhmf2(doy,latitude,height,elev);

%     Routine to compute the hydrostatic mapping function nhmf2 which
%     depends on DOY (day of year) and station position (latitude 
%     and height above geoid).
 
% Define variables used
%   a,b,c       - the a,b,and c coeffiecents in the continued fraction
%                 form of Marini
%   beta        - intermediate term in calculation
%   gamma       - intermediate term in calculation
%   sine        - sine of elevation angle
%   cose        - cos of elevation angle
%   hmf(1)      - delay mapping function
%   hmf(2)      - d_mapping_function/d_elevation (dhmf2/d_el)
%   topcon      - constant of top of mapping fuinction to ensure
%                 that value is 1.0000 at zenith 
%   height     - height of site above geoid (meters)
%   hs_km      - Height of site in kms.
%   latitude   - latitude (degrees)
%   l          - absolute latitude
%   dl         - incremental latitude from last lat_hmf
%   elev       - elevation (degrees)
%   epoch      - if Julian date of observation is known for the observation,
%              - then epoch can be used to get day of year.
%              - (if epoch  is passed as argument, then un-comment the
%                 line converting epoch to doy.)
%   doy        - days since Dec 31 
%   doy_atm    - doy for atmosphere relative to Jan 28.
%   doyr_atm   - doy_atm in radians;
%   cost       - cosine(day of year)
%   doy2rad    - convert doy to radians
%   lat_hmf     - latitudes at which coefficients are defined (5).
%   abc_avg     - continued fraction coefficients at latitudes lat_hmf
%   abc_amp     - amplitude of annual variation of abc_avg
%   daavg, daamp, etc - incremental values for interpolation
%   aavg,  aamp,  etc - average and amplitude at latitude
%   a_ht, b_ht, c_ht - parameters for continued fraction for height corr'n.
%   define parameters used for calculating coefficients.

lat_hmf = [15, 30, 45, 60, 75];

abc_avg = [ ...
     1.2769934e-3,1.2683230e-3,1.2465397e-3,1.2196049e-3,1.2045996e-3; ...
     2.9153695e-3,2.9152299e-3,2.9288445e-3,2.9022565e-3,2.9024912e-3; ...
     62.610505e-3,62.837393e-3,63.721774e-3,63.824265e-3,64.258455e-3]';

abc_amp = [ ... 
     0.0,   1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5; ...
     0.0,   2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5; ...
     0.0,   9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5]';

a_ht = 2.53e-5;
b_ht = 5.49e-3;
c_ht = 1.14e-3;

%   conversions:
doy2rad = 2*3.14159265/365.25;
deg2rad = 3.14159265/180.0;

%   convert height in meters to kilometers
hs_km  = height/1000.0;

l = abs(latitude);
if (latitude < 0) 
   doy = doy + 365.25/2;
end

% mod aen 930517 Use phase of 28 days (winter extremum corresponds to Jan 28)
%                based on least-square fit to 
%                raytrace of radiosonde data for DRT, ELP, ALB, CHH, FAI,
%                MUN, and LIH.
doy_atm  = doy - 28.0;
doyr_atm = doy_atm * doy2rad;
cost = cos(doyr_atm);

%   Coefficients for the continued fraction expansion for each latitude.
%   for latitudes less than 15 degrees:
I_low = find(l <= lat_hmf(1));
if ~isempty(I_low)
	a(I_low) = abc_avg(1,1);
	b(I_low) = abc_avg(1,2);
	c(I_low) = abc_avg(1,3);
end % if

%   for latitudes between 15 and 75  degrees:
for i = 1:4
	I_mid = find(l > lat_hmf(i) & l <= lat_hmf(i+1));
   if ~isempty(I_mid)
	   dl = (l-lat_hmf(i))/(lat_hmf(i+1)-lat_hmf(i));
   	daavg =   abc_avg(i+1,1)-abc_avg(i,1);
		daamp =   abc_amp(i+1,1)-abc_amp(i,1);
		aavg  =   abc_avg(i,1) + dl(I_mid)*daavg;
		aamp  =   abc_amp(i,1) + dl(I_mid)*daamp;
      a(I_mid) = aavg - aamp.*cost(I_mid);

		dbavg =   abc_avg(i+1,2)-abc_avg(i,2);
		dbamp =   abc_amp(i+1,2)-abc_amp(i,2);
		bavg  =   abc_avg(i,2) + dl(I_mid)*dbavg;
		bamp  =   abc_amp(i,2) + dl(I_mid)*dbamp;
		b(I_mid) = bavg - bamp.*cost(I_mid);

		dcavg =   abc_avg(i+1,3)-abc_avg(i,3);
		dcamp =   abc_amp(i+1,3)-abc_amp(i,3);
		cavg  =   abc_avg(i,3) + dl(I_mid)*dcavg;
		camp  =   abc_amp(i,3) + dl(I_mid)*dcamp; 
		c(I_mid)     = cavg - camp.*cost(I_mid);
	end % if
end % for

%   for latitudes greater than 75 degrees:
I_high = find(l > lat_hmf(5));
if ~isempty(I_high)
	a(I_high) = abc_avg(5,1) - abc_amp(5,1)*cost;
	b(I_high) = abc_avg(5,2) - abc_amp(5,2)*cost;
	c(I_high) = abc_avg(5,3) - abc_amp(5,3)*cost;
end % if

%   Now the coefficients exist; calculate the mapping function, hmf(1),
%       and the change of mapping function with elevation, 
%       dhmf/d_el = hmf(2).
%   To get delay-rate correction d_tau/dt:
%      d_tau/dt = d_tau-zen/dt*hmf(1) + tau-zen*dhmf/d_el*d_el/dt
sine   = sin(elev * deg2rad);
cose   = cos(elev * deg2rad);
beta   = b' ./ (sine + c' );
gamma  = (a' ./ (sine + beta));
topcon = (1.0 + a ./ (1.0 + b ./ (1.0 + c)))';

hmf(:,1) = topcon ./ (sine + gamma);

hmf(:,2) = -topcon ./ (sine + gamma).^2 .* ...
	(cose - a' ./ (sine + beta).^2 .* cose .* ...
	(1.0 - b' ./ (sine + c').^2 ));

%   Apply height correction to mapping function only 
%         (not to dmf/d_el since this is a small correction): 
%      1) height correction coefficient is 
%         1/sine(elev) - continued fraction(a_ht,b_ht,c_ht).
%      2) height correction is ht_corr_coef times height in km.
beta   = b_ht ./ (sine + c_ht);
gamma  = a_ht ./ ( sine + beta);
topcon = (1.0 + a_ht ./ (1.0 + b_ht ./ (1.0 + c_ht)));
ht_corr_coef = 1 ./ sine - topcon ./ (sine + gamma);
ht_corr = ht_corr_coef; 					% hs_km
hmf(:,1) = hmf(:,1) + ht_corr;

% end of nhmf2

%%%%% END ALGORITHM CODE %%%%%

% end of TROPUNB1
