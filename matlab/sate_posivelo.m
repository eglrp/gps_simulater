function [satp,satv]=sate_posivelo(t,orbit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Scope:  �ɹ���������������ECEF�����е�λ��
% Usage:  P=sate_position(t,orbit)
% Description of parameter:
%      t               -input,t ���յ��õ���ʱ�Ľ��ջ�GPSʱ��
%      orbit(1)    -input,Crs ���ǵ��ľ�ص��͸����������������
%      orbit(2)    -input,delt_n  ƽ�������ٶȲ�
%      orbit(3)    -input,M0  �ο�ʱ�̵�ƽ�����
%      orbit(4)    -input,Cuc ������ǵĵ��͸��������,������
%      orbit(5)    -input,e ��Բ���������
%      orbit(6)    -input,Cus ������ǵĵ��͸��������,������
%      orbit(7)    -input,sqrta ���������ƽ����
%      orbit(8)    -input,toe �ο���Ԫ
%
%      orbit(9)    -input,Cic �����ǵĵ��͸����������������
%      orbit(10)   -input,omega0 ������׼�ྭ
%      orbit(11)   -input,Cis �����ǵĵ��͸����������������
%      orbit(12)   -input,i0 ������
%      orbit(13)   -input,Crc ���ǵ��ľ�ص��͸����������������
%      orbit(14)   -input,w0  �ο�ʱ�̵Ľ��ص�Ǿ�
%      orbit(15)   -input,delt_omega ������ྭ�仯��
%      orbit(16)   -input,d_i    �����Ǳ仯��di/dt

%      orbit(17)   -input,af0     ����ʱ��ƫ��
%      orbit(18)   -input,af1     ����ʱ��Ư��
%      orbit(19)   -input,af2     ����ʱ��Ư����
%      orbit(20)   -input,toc     ʱ��ʱ��
%      pole        -input,      �����Ʋ���
%      P           -output,P  ������ECEF����ϵ�µ�λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% constant��
GM    = 3.986005e14;  %������������
wie   = 7.2921151467e-5;  %������ת���ٶ�
      
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
% ��������ƽ��������
n=sqrt(GM/a^3)+delt_n;
% ����黯ʱ��delt_t
delt_t=t-toe;
% ����۲���Ԫt��ƽ�����M
M=M0+n*delt_t;
% ����ƫ�����E
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
% �������ǵĵ���ʸ��r
r=a*(1-e*cos(E));
% ����������f
f=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
% ����������Ǿ�fai
fai0=w0+f;
% �����㶯������
delt_u=Cus*sin(2*fai0)+Cuc*cos(2*fai0);
delt_r=Crs*sin(2*fai0)+Crc*cos(2*fai0);
delt_i=Cis*sin(2*fai0)+Cic*cos(2*fai0);
% ���㾭���㶯������������Ǿ࣬����ʸ���͹�������
fai=fai0+delt_u;
r=r+delt_r;
i=i0+delt_i+d_i*delt_t;
%����۲���Ԫt�������㾭��na
na=omega0+(delt_omega-wie)*delt_t-wie*toe;
% ���������ڹ��ƽ��ֱ������ϵ�е�λ��
X0=r*[cos(fai) sin(fai) 0]';
% ����������Э���������ϵ��λ��
C=[cos(na) -sin(na)*cos(i) sin(na)*sin(i);sin(na) cos(na)*cos(i) -cos(na)*sin(i);0 sin(i) cos(i)];
satp=C*X0;
% Cpole=[1 0 xp;0 1 yp;-xp,-yp 1];
% P=Cpole*P;

%%
% ƫ�������d_E
d_E=n/(1-e*cos(E));
% ������Ǿ��󵼣�У��ǰ��d_fai0
d_fai0=d_E*sqrt((1+e)/(1-e))*(cos(f/2))^2/(cos(E/2))^2;
% �㶯��������
d_u=2*d_fai0*(Cus*cos(2*fai0)-Cuc*sin(2*fai0));
d_r=2*d_fai0*(Crs*cos(2*fai0)-Crc*sin(2*fai0));
d_ii=2*d_fai0*(Cis*cos(2*fai0)-Cic*sin(2*fai0));
% ������ǾࣨУ�����󵼣�����ʸ���󵼣���������
d_fai=d_fai0+d_u;
d_r=a*e*sin(E)*d_E+d_r;
d_i=delt_i+d_ii;
% ����۲���Ԫt�������㾭�ȵ���d_na
d_na=delt_omega-wie;
% ���������ڹ��ֱ������ϵ�е��ٶ�
v0=d_r*[cos(fai) sin(fai) 0]'+r*d_fai*[-sin(fai) cos(fai) 0]';
% ����������˲ʱ����ϵ���ٶ�
C2=[-d_na*sin(na)  -d_na*cos(na)*cos(i)+d_i*sin(na)*sin(i)   d_na*cos(na)*sin(i)+d_i*sin(na)*cos(i);...
     d_na*cos(na)  -d_na*sin(na)*cos(i)-d_i*cos(na)*sin(i)   d_na*sin(na)*sin(i)-d_i*cos(na)*cos(i);...
     0              d_i*cos(i)                              -d_i*sin(i)];
satv=C*v0+C2*X0;
% ����������Э�����ϵ�е��ٶ�
% % % V=Cpole*V;



