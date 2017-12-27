%
%
%    地理系速度向地球系系速度转换    Navi——>ECEF
%
%
%
%=========================================================================

function  veloE = veloEN(posi,veloN)

long=posi(1,1)*pi/180.0;
lati=posi(2,1)*pi/180.0;

    % 坐标系转换矩阵E-->N(东北天）
Cne=[-sin(long),          cos(long),           0;
     -sin(lati)*cos(long),  -sin(lati)*sin(long),    cos(lati);
      cos(lati)*cos(long),   cos(lati)*sin(long),    sin(lati)];
  
veloE=Cne'*veloN;