%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    用户位置速度计算
%
%  输入参数:
%  t          仿真时间
%  T          仿真步长
%  atti       横滚、俯仰、航向（单位：度）
%  atti_rate  横滚速率、俯仰速率、航向速率（单位：度/秒）
%  veloB      飞机运动速度——X右翼、Y机头、Z天向（单位：米/秒）
%  acceB      飞机运动加速度——X右翼、Y机头、Z天向（单位：米/秒/秒）
%  posi       航迹发生器初始位置经度、纬度、高度（单位：度、度、米）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xr,Vr]=user_positions(t,T,atti,atti_rate,veloB,acceB,posi)

Re=6378137.0;                   %地球半径     （单位：米） 
f=1/298.257223563;          %地球的椭圆扁率

    %飞行器位置
long=posi(1,1)*pi/180.0;lati=posi(2,1)*pi/180.0;heig=posi(3,1);    
    %地球曲率半径求解
Rm=Re*(1-2*f+3*f*sin(lati)*sin(lati));
Rn=Re*(1+f*sin(lati)*sin(lati));
    
    % 航迹仿真
 %[t,atti,atti_rate,veloB,acceB]=trace_velo(t,T,atti,atti_rate,veloB,acceB);

    % 姿态角
roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
    
    % 坐标系转换矩阵N-->B
Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
     cos(pitch)*sin(head),                                cos(pitch)*cos(head),                                sin(pitch);
     sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
    
    % 地理系速度、加速度
veloN=Cbn'*veloB;
acceN=Cbn'*acceB;
Ve=veloN(1,1);Vn=veloN(2,1);Vu=veloN(3,1);
    
    % 地理系位置（经纬高）
heig=heig+t*Vu;
lati=lati+t*(Vn/(Rm+heig));
long=long+t*(Ve/((Rn+heig)*cos(lati)));
    
posi(1,1)=long*180.0/pi;
posi(2,1)=lati*180.0/pi;
posi(3,1)=heig;
posiE = geo2ecef(posi);  

    % 坐标系转换矩阵ECEF-->N(东北天）
Cne=[-sin(long),          cos(long),           0;
     -sin(lati)*cos(long),  -sin(lati)*sin(long),    cos(lati);
      cos(lati)*cos(long),   cos(lati)*sin(long),    sin(lati)];
      
    % 地球系速度、位置（经纬高）、加速度
veloE=Cne'*veloN;

Xr=posiE;
Vr=veloE;
% acceE=Cne'*acceN;
    

    