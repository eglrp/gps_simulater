function [XS_tx, VS_tx, time_M, time_R] = satellite_positions(t,orbit,nsat)
%	compute satellite positions and transit time
%	Input:  
%           t                       �źŽ���ʱ��
%		    orbit                 ���ǹ����Ϣ
%	Output:  
%	        XS_tx   VS_tx    �źŷ���ʱ��������ECEF�µ�λ�á��ٶ� 
%           time_M           ���ǵ�����M�Ĵ���ʱ��
%           time_R            ���ǵ�����R�Ĵ���ʱ��

%% Initialize ======================================================
global sign_set;

time_M=zeros(nsat,1);
XS_tx=zeros(nsat,3);
VS_tx=zeros(nsat,3);
% position in GEOREF : long (deg)  lati (deg)  alti (m)
posiM		= sign_set.init_posi_M;
posiR		= sign_set.init_posi_R;

%% Master 
traceM.posi_GEO	= posiM;
traceM.posi_ECEF = geo2ecef(posiM);

%% Rover 
traceR.posi_GEO = posiR;
traceR.posi_ECEF = geo2ecef(posiR);

%% �����վM������ʱ
Xr1 = traceM.posi_ECEF;
t_orbit = sign_set.t_sate+t;        
for j=1:nsat
    % ��������λ���ٶ�
     sprintf('%d',t_orbit),
     sprintf('%d',orbit(j,:)),
    [Xs,Vs]=sate_posivelo(t_orbit,orbit(j,:));      %�źŽ���ʱ������λ��
     sprintf('%d',Xs),
     sprintf('%d',Vs),

    %% �������㴫����ʱ
    Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
    Tp_old1=0;
    while((Tp1-Tp_old1)>1e-12)
        % ����У��������λ��
        [Xs,Vs]=sate_posivelo(t_orbit-Tp1,orbit(j,:));            %�źŷ���ʱ������λ��
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
    [XS_tx(j,:),VS_tx(j,:)]=sate_posivelo(t_orbit-time_M(j,1),orbit(j,:));   
end

% �������վR�Ĵ�����ʱ
Xr2 = traceR.posi_ECEF;
for j=1:nsat
    % ��������λ���ٶ�
    [Xs,Vs]=sate_posivelo(t_orbit,orbit(j,:));
    %% �������㴫����ʱ
    Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
    Tp_old2=0;
    while((Tp2-Tp_old2)>1e-12)
        % ����У��������λ��
        [Xs,Vs]=sate_posivelo(t_orbit-Tp2,orbit(j,:));
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
    [XS_tx(j,:),VS_tx(j,:)]=sate_posivelo(t_orbit-time_R(j,1),orbit(j,:));   
end

%XS_tx, VS_tx, time_M, time_R

