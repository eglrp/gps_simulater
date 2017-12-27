%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    �û�λ���ٶȼ���
%
%  �������:
%  t          ����ʱ��
%  T          ���沽��
%  atti       ��������������򣨵�λ���ȣ�
%  atti_rate  ������ʡ��������ʡ��������ʣ���λ����/�룩
%  veloB      �ɻ��˶��ٶȡ���X����Y��ͷ��Z���򣨵�λ����/�룩
%  acceB      �ɻ��˶����ٶȡ���X����Y��ͷ��Z���򣨵�λ����/��/�룩
%  posi       ������������ʼλ�þ��ȡ�γ�ȡ��߶ȣ���λ���ȡ��ȡ��ף�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xr,Vr]=user_positions(t,T,atti,atti_rate,veloB,acceB,posi)

Re=6378137.0;                   %����뾶     ����λ���ף� 
f=1/298.257223563;          %�������Բ����

    %������λ��
long=posi(1,1)*pi/180.0;lati=posi(2,1)*pi/180.0;heig=posi(3,1);    
    %�������ʰ뾶���
Rm=Re*(1-2*f+3*f*sin(lati)*sin(lati));
Rn=Re*(1+f*sin(lati)*sin(lati));
    
    % ��������
 %[t,atti,atti_rate,veloB,acceB]=trace_velo(t,T,atti,atti_rate,veloB,acceB);

    % ��̬��
roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
    
    % ����ϵת������N-->B
Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
     cos(pitch)*sin(head),                                cos(pitch)*cos(head),                                sin(pitch);
     sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
    
    % ����ϵ�ٶȡ����ٶ�
veloN=Cbn'*veloB;
acceN=Cbn'*acceB;
Ve=veloN(1,1);Vn=veloN(2,1);Vu=veloN(3,1);
    
    % ����ϵλ�ã���γ�ߣ�
heig=heig+t*Vu;
lati=lati+t*(Vn/(Rm+heig));
long=long+t*(Ve/((Rn+heig)*cos(lati)));
    
posi(1,1)=long*180.0/pi;
posi(2,1)=lati*180.0/pi;
posi(3,1)=heig;
posiE = geo2ecef(posi);  

    % ����ϵת������ECEF-->N(�����죩
Cne=[-sin(long),          cos(long),           0;
     -sin(lati)*cos(long),  -sin(lati)*sin(long),    cos(lati);
      cos(lati)*cos(long),   cos(lati)*sin(long),    sin(lati)];
      
    % ����ϵ�ٶȡ�λ�ã���γ�ߣ������ٶ�
veloE=Cne'*veloN;

Xr=posiE;
Vr=veloE;
% acceE=Cne'*acceN;
    

    