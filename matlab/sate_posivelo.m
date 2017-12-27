function [satp,satv]=sate_posivelo(t,orbit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Scope:  由轨道参数求解卫星在ECEF坐标中的位置
% Usage:  P=sate_position(t,orbit)
% Description of parameter:
%      t               -input,t 接收到该电文时的接收机GPS时间
%      orbit(1)    -input,Crs 卫星地心距地调和改正项振幅，正弦项
%      orbit(2)    -input,delt_n  平均运行速度差
%      orbit(3)    -input,M0  参考时刻的平近点角
%      orbit(4)    -input,Cuc 升交距角的调和改正项振幅,余弦项
%      orbit(5)    -input,e 椭圆轨道离心率
%      orbit(6)    -input,Cus 升交距角的调和改正项振幅,正弦项
%      orbit(7)    -input,sqrta 轨道长半轴平方根
%      orbit(8)    -input,toe 参考历元
%
%      orbit(9)    -input,Cic 轨道倾角的调和改正项振幅，余弦项
%      orbit(10)   -input,omega0 升交点准赤经
%      orbit(11)   -input,Cis 轨道倾角的调和改正项振幅，正弦项
%      orbit(12)   -input,i0 轨道倾角
%      orbit(13)   -input,Crc 卫星地心距地调和改正项振幅，余弦项
%      orbit(14)   -input,w0  参考时刻的近地点角距
%      orbit(15)   -input,delt_omega 升交点赤经变化率
%      orbit(16)   -input,d_i    轨道倾角变化率di/dt

%      orbit(17)   -input,af0     卫星时钟偏差
%      orbit(18)   -input,af1     卫星时钟漂移
%      orbit(19)   -input,af2     卫星时钟漂移率
%      orbit(20)   -input,toc     时钟时间
%      pole        -input,      地球极移参数
%      P           -output,P  卫星在ECEF坐标系下的位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% constant：
GM    = 3.986005e14;  %地球引力常数
wie   = 7.2921151467e-5;  %地球自转角速度
      
Crs         = orbit(1);
delt_n      = orbit(2);
M0          = orbit(3);
Cuc         = orbit(4);
e           = orbit(5);
Cus         = orbit(6);
a           = orbit(7)^2;
toe         = orbit(8);

Cic         = orbit(9);
omega0      = orbit(10);
Cis         = orbit(11);
i0          = orbit(12);

Crc        = orbit(13);
w0         = orbit(14);
delt_omega = orbit(15);
d_i        = orbit(16);
af0        = orbit(17);
af1        = orbit(18);
af2        = orbit(19);
toc        = orbit(20);

% xp         = pole(1);
% yp         = pole(2);
% 卫星运行平均角速率
n=sqrt(GM/a^3)+delt_n;
% 计算归化时间delt_t
delt_t=t-toe;
% 计算观测历元t的平近点角M
M=M0+n*delt_t;
% 计算偏近点角E
E0=M;
E=E0;
delt_M=1;
while  ~(abs(delt_M)<(1e-12))
    M1=E0-e*sin(E0);
    delt_M=M-M1;
    delt_E=delt_M/(1-e*cos(E0));
    E=E0+delt_E;
    E0=E;
end
% 计算卫星的地心矢径r
r=a*(1-e*cos(E));
% 计算真近点角f
f=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
% 计算升交点角距fai
fai0=w0+f;
% 计算摄动改正项
delt_u=Cus*sin(2*fai0)+Cuc*cos(2*fai0);
delt_r=Crs*sin(2*fai0)+Crc*cos(2*fai0);
delt_i=Cis*sin(2*fai0)+Cic*cos(2*fai0);
% 计算经过摄动改正的升交点角距，卫星矢径和轨道面倾角
fai=fai0+delt_u;
r=r+delt_r;
i=i0+delt_i+d_i*delt_t;
%计算观测历元t的升交点经度na
na=omega0+(delt_omega-wie)*delt_t-wie*toe;
% 计算卫星在轨道平面直角坐标系中的位置
X0=r*[cos(fai) sin(fai) 0]';
% 计算卫星在协议地球坐标系的位置
C=[cos(na) -sin(na)*cos(i) sin(na)*sin(i);sin(na) cos(na)*cos(i) -cos(na)*sin(i);0 sin(i) cos(i)];
satp=C*X0;
% Cpole=[1 0 xp;0 1 yp;-xp,-yp 1];
% P=Cpole*P;

%%
% 偏近点角求导d_E
d_E=n/(1-e*cos(E));
% 升交点角距求导（校正前）d_fai0
d_fai0=d_E*sqrt((1+e)/(1-e))*(cos(f/2))^2/(cos(E/2))^2;
% 摄动改正项求导
d_u=2*d_fai0*(Cus*cos(2*fai0)-Cuc*sin(2*fai0));
d_r=2*d_fai0*(Crs*cos(2*fai0)-Crc*sin(2*fai0));
d_ii=2*d_fai0*(Cis*cos(2*fai0)-Cic*sin(2*fai0));
% 升交点角距（校正后）求导，地心矢径求导，轨道倾角求导
d_fai=d_fai0+d_u;
d_r=a*e*sin(E)*d_E+d_r;
d_i=delt_i+d_ii;
% 计算观测历元t的升交点经度导数d_na
d_na=delt_omega-wie;
% 计算卫星在轨道直角坐标系中的速度
v0=d_r*[cos(fai) sin(fai) 0]'+r*d_fai*[-sin(fai) cos(fai) 0]';
% 计算卫星在瞬时地球系的速度
C2=[-d_na*sin(na)  -d_na*cos(na)*cos(i)+d_i*sin(na)*sin(i)   d_na*cos(na)*sin(i)+d_i*sin(na)*cos(i);...
     d_na*cos(na)  -d_na*sin(na)*cos(i)-d_i*cos(na)*sin(i)   d_na*sin(na)*sin(i)-d_i*cos(na)*cos(i);...
     0              d_i*cos(i)                              -d_i*sin(i)];
satv=C*v0+C2*X0;
% 计算卫星在协议地球系中的速度
% % % V=Cpole*V;



