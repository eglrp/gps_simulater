function [time_M, time_R, ex_M,ey_M,ez_M,ex_R, ey_R, ez_R] = TimeDelay_SR(orbit, Cnsat, Xr_M, Xr_R)
%	compute satellite positions and transit time
%	Input:  
%           t                        信号接收时间
%		    orbit                 卫星轨道信息
%           Xr_M, Xr_R       载体位置
%           Cnsat                可见卫星
%	Output:  
%           time_M           卫星到载体M的传输时间
%           time_R            卫星到载体R的传输时间
%          ex_M, ey_M, ez_M        卫星到载体M方向余弦
%          ex_R, ey_R,  ez_R           卫星到载体R方向余弦
%% Initialize ======================================================
global sign_set;

for t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time
    nsat0=size(Cnsat,1);                 % 可见星数目
    i=round(1+5*t);                        % 采样数据点数
%% Master 
traceM.posi_ECEF = Xr_M;
%% Rover 
traceR.posi_ECEF = Xr_R;
%% 求解主站M传输延时
Xr1 = traceM.posi_ECEF;
t_orbit = sign_set.t_sate+t;    
        for j=1:nsat0
            % 解算卫星位置速度
           % j,
            k=Cnsat(j,1);
            [Xs,Vs]=sate_posivelo(t_orbit,orbit(k,:));      %信号接收时刻卫星位置
            %% 迭代计算传输延时
            Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
            Tp_old1=0;
            while((Tp1-Tp_old1)>1e-12)
                % 计算校正后卫星位置
                [Xs,Vs]=sate_posivelo(t_orbit-Tp1,orbit(k,:));            %信号发射时刻卫星位置
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
            time_M(j,i)=Tp1;
            % 求解卫星到主站的单位方向余弦
            ex_M(j,i)= (Xs(1)-Xr1(1))/sqrt((Xs-Xr1)'*(Xs-Xr1));      %  x坐标
            ey_M(j,i)= (Xs(2)-Xr1(2))/sqrt((Xs-Xr1)'*(Xs-Xr1));      %  y坐标
            ez_M(j,i)= (Xs(3)-Xr1(3))/sqrt((Xs-Xr1)'*(Xs-Xr1));      %  z坐标  
        end
        % 求解流动站R的传输延时
        Xr2 = traceR.posi_ECEF;
        for j=1:nsat0
            % j,
            k=Cnsat(j,1);
            % 解算卫星位置速度
            [Xs,Vs]=sate_posivelo(t_orbit,orbit(k,:));
            %% 迭代计算传输延时
            Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
            Tp_old2=0;
            while((Tp2-Tp_old2)>1e-12)
                % 计算校正后卫星位置
                [Xs,Vs]=sate_posivelo(t_orbit-Tp2,orbit(k,:));
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
            time_R(j,i)=Tp2;
            % 求解卫星到流动站的单位方向余弦
            ex_R(j,i)= (Xs(1)-Xr2(1))/sqrt((Xs-Xr2)'*(Xs-Xr2));      %  x坐标
            ey_R(j,i)= (Xs(2)-Xr2(2))/sqrt((Xs-Xr2)'*(Xs-Xr2));      %  y坐标
            ez_R(j,i)= (Xs(3)-Xr2(3))/sqrt((Xs-Xr2)'*(Xs-Xr2));      %  z坐标  
        end
end
