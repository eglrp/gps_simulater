function posi_ECEF = geo2ecef(posi_GEO)
%GEO2CART Conversion of geographical coordinates (lambda,phi,  h) to
%Cartesian coordinates (X, Y, Z). 
%[X, Y, Z] = geo2cart(lambda,phi,  h);
%Format for phi and lambda: [degrees minutes seconds].
%h, X, Y, and Z are in meters.
%   Inputs:       
%       lambda    - geocentric longitude (format [degrees minutes seconds]) 
%       phi           - geocentric latitude (format [degrees minutes seconds])
%       h              - height
%   Outputs:
%       X, Y, Z   - Cartesian coordinates (meters)
%==========================================================================
% constant

Re =  6378137;                          % GPS (WGS-84)      Ellipsoid semi-major axis [m]
f =1/298.257223563;                  % GPS (WGS-84)      Ellipsoid flattening
Wie  = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
g = 9.7803698;

% user position (rad)
long = posi_GEO(1,1) * pi / 180.0;
lati = posi_GEO(2,1) * pi / 180.0;
heig = posi_GEO(3,1);

% earth ellipticity
Rn = Re * (1 + f * sin(lati) * sin(lati));

% position ECEF
X = (Rn + heig) * cos(lati) * cos(long); 
Y = (Rn + heig) * cos(lati) * sin(long) ;
Z = (Rn * (1 - f)^2 + heig) * sin(lati);
posi_ECEF = [X Y Z]';

%%%%%%%%%%%%%% end geo2cart.m  %%%%%%%%%%%%%%%%%%%%%%%%
