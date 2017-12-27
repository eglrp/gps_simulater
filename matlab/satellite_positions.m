function [XS_tx, VS_tx, time_M, time_R] = satellite_positions(t,orbit,nsat)
%	compute satellite positions and transit time
%	Input:  
%           t                       信号接收时间
%		    orbit                 卫星轨道信息
%	Output:  
%	        XS_tx   VS_tx    信号发射时刻卫星在ECEF下的位置、速度 
%           time_M           卫星到载体M的传输时间
%           time_R            卫星到载体R的传输时间

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

%% 求解主站M传输延时
Xr1 = traceM.posi_ECEF;
t_orbit = sign_set.t_sate+t;        
for j=1:nsat
    % 解算卫星位置速度
     sprintf('%d',t_orbit),
     sprintf('%d',orbit(j,:)),
    [Xs,Vs]=sate_posivelo(t_orbit,orbit(j,:));      %信号接收时刻卫星位置
     sprintf('%d',Xs),
     sprintf('%d',Vs),

    %% 迭代计算传输延时
    Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
    Tp_old1=0;
    while((Tp1-Tp_old1)>1e-12)
        % 计算校正后卫星位置
        [Xs,Vs]=sate_posivelo(t_orbit-Tp1,orbit(j,:));            %信号发射时刻卫星位置
        %地球旋转改正
        % 由t_orbit-Tp时刻的ECEF坐标系转换到t_orbit时刻的ECEF坐标系     
        C=[ cos(sign_set.wie*Tp1)    sin(sign_set.wie*Tp1)    0;...
            -sin(sign_set.wie*Tp1)    cos(sign_set.wie*Tp1)    0;...
            0              0              1];
        Xs=C*Xs;
        Vs=C*Vs;
        % 修正传输延时
        Tp_old1=Tp1;
        Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
    end
    time_M(j,1)=Tp1;
    [XS_tx(j,:),VS_tx(j,:)]=sate_posivelo(t_orbit-time_M(j,1),orbit(j,:));   
end

% 求解流动站R的传输延时
Xr2 = traceR.posi_ECEF;
for j=1:nsat
    % 解算卫星位置速度
    [Xs,Vs]=sate_posivelo(t_orbit,orbit(j,:));
    %% 迭代计算传输延时
    Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
    Tp_old2=0;
    while((Tp2-Tp_old2)>1e-12)
        % 计算校正后卫星位置
        [Xs,Vs]=sate_posivelo(t_orbit-Tp2,orbit(j,:));
        % 由t_orbit-Tp时刻的ECEF坐标系转换到t_orbit时刻的ECEF坐标系
        C=[ cos(sign_set.wie*Tp2)    sin(sign_set.wie*Tp2)    0;...
            -sin(sign_set.wie*Tp2)    cos(sign_set.wie*Tp2)    0;...
            0              0              1];
        Xs=C*Xs;
        Vs=C*Vs;
        % 修正传输延时
        Tp_old2=Tp2;
        Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
    end
    time_R(j,1)=Tp2;
    [XS_tx(j,:),VS_tx(j,:)]=sate_posivelo(t_orbit-time_R(j,1),orbit(j,:));   
end

%XS_tx, VS_tx, time_M, time_R

