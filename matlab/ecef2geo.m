function posi_GEO = ecef2geo(posi_ECEF)
%CART2GEO Conversion of Cartesian coordinates (X,Y,Z) to geographical
%coordinates (lambda,phi,  h) on a selected reference ellipsoid.
%将直角坐标转换为经纬高坐标
%
%==========================================================================

% constant
Re =  6378137;                          % GPS (WGS-84)      Ellipsoid semi-major axis [m]
f =1/298.257223563;                  % GPS (WGS-84)      Ellipsoid flattening
Wie  = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
g = 9.7803698;


% user position
X = posi_ECEF(1);
Y = posi_ECEF(2);
Z = posi_ECEF(3);

e = sqrt(1-(1-f)^2);

% position GEO
long = atan2(Y,X);   %经度
heig = 0.1; oldh = 0;               %高度初值设为0
lati = atan(Z/((1-f^2)*sqrt(X^2+Y^2)));  %纬度初值 高度初值设为0
i=0;
while abs(heig-oldh) > 1.e-10
	oldh = heig;
	RpH = X/cos(lati)/cos(long);   %lati 
	Rn = Re*(1 + f*sin(lati)^2);
% 	Rn = Re/sqrt(cos(lati)^2+(1-e*e)*sin(lati)^2);
	lati = atan(RpH/(RpH-Rn*e*e)*Z/sqrt(X^2+Y^2));  %经度
	heig = RpH - Rn;    
	i = i+1;
	if i>200
		disp('ECEF -> GEO Error!');
		break;
	end
end

posi_GEO = [long*180/pi lati*180/pi heig]';

%%%%%%%%%%%%%% end cart2geo.m %%%%%%%%%%%%%%%%%%%
