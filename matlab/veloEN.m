%
%
%    ����ϵ�ٶ������ϵϵ�ٶ�ת��    Navi����>ECEF
%
%
%
%=========================================================================

function  veloE = veloEN(posi,veloN)

long=posi(1,1)*pi/180.0;
lati=posi(2,1)*pi/180.0;

    % ����ϵת������E-->N(�����죩
Cne=[-sin(long),          cos(long),           0;
     -sin(lati)*cos(long),  -sin(lati)*sin(long),    cos(lati);
      cos(lati)*cos(long),   cos(lati)*sin(long),    sin(lati)];
  
veloE=Cne'*veloN;