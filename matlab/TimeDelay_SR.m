function [time_M, time_R] = TimeDelay_SR(t, orbit, Cnsat, Xr_M, Xr_R)
%	compute satellite positions and transit time
%	Input:  
%           t                        �źŽ���ʱ��
%		    orbit                 ���ǹ����Ϣ
%           Xr_M, Xr_R       ����λ��
%           Cnsat                �ɼ�����
%	Output:  
%           time_M           ���ǵ�����M�Ĵ���ʱ��
%           time_R            ���ǵ�����R�Ĵ���ʱ��

%% Initialize ======================================================
global sign_set;
nsat0=size(Cnsat,1);                 % �ɼ�����Ŀ
time_M=zeros(nsat0,1);
time_R=zeros(nsat0,1);
%% Master 
traceM.posi_ECEF = Xr_M;
%% Rover 
traceR.posi_ECEF = Xr_R;
%% �����վM������ʱ
Xr1 = traceM.posi_ECEF;
t_orbit = sign_set.t_sate+t;    

for j=1:nsat0
    % ��������λ���ٶ�
   % j,
    i=Cnsat(j,1);
    [Xs,Vs]=sate_posivelo(t_orbit,orbit(i,:));      %�źŽ���ʱ������λ��
    %% �������㴫����ʱ
    Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
    Tp_old1=0;
    while((Tp1-Tp_old1)>1e-12)
        % ����У��������λ��
        [Xs,Vs]=sate_posivelo(t_orbit-Tp1,orbit(i,:));            %�źŷ���ʱ������λ��
        %������ת����
        % ��t_orbit-Tpʱ�̵�ECEF����ϵת����t_orbitʱ�̵�ECEF����ϵ     
        C=[ cos(sign_set.wie*Tp1)    sin(sign_set.wie*Tp1)    0;...
            -sin(sign_set.wie*Tp1)    cos(sign_set.wie*Tp1)    0;...
            0              0              1];
        Xs=C*Xs;
        Vs=C*Vs;
        % ����������ʱ
        Tp_old1=Tp1;
        Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
    end
    time_M(j,1)=Tp1;
end

% �������վR�Ĵ�����ʱ
Xr2 = traceR.posi_ECEF;
for j=1:nsat0
    % j,
    i=Cnsat(j,1);
    % ��������λ���ٶ�
    [Xs,Vs]=sate_posivelo(t_orbit,orbit(i,:));
    %% �������㴫����ʱ
    Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
    Tp_old2=0;
    while((Tp2-Tp_old2)>1e-12)
        % ����У��������λ��
        [Xs,Vs]=sate_posivelo(t_orbit-Tp2,orbit(i,:));
        % ��t_orbit-Tpʱ�̵�ECEF����ϵת����t_orbitʱ�̵�ECEF����ϵ
        C=[ cos(sign_set.wie*Tp2)    sin(sign_set.wie*Tp2)    0;...
            -sin(sign_set.wie*Tp2)    cos(sign_set.wie*Tp2)    0;...
            0              0              1];
        Xs=C*Xs;
        Vs=C*Vs;
        % ����������ʱ
        Tp_old2=Tp2;
        Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
    end
    time_R(j,1)=Tp2;
end


