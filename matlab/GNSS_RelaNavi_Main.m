 %	Main funciton of GNSS_RelaNavi_Main
%	Input:
%		void
%	Output:
%		void
%
%	Created by Wu Ling  from NUAA NRC at 2017/10/06
%==========================================================================
%% %%清除所有变量、关闭所有打开的图%%%%%%%%
clear all;
close all;
clc;

%% %%%%%%%%%Get initial settings %%%%%%%%%%%%%%%%%
global sign_set
%  获得初始化参数
sign_set = Inital_Settings();     
lamada1=sign_set.GPS_bo1;
lamada2=sign_set.GPS_bo2;
Deg0=sign_set.cutoff;
%posi_GEO    =zeros(3,sign_set.datalength);      %地理坐标系下位置
T = sign_set.TDperiod;% 航迹生成(传输延时）周期
% Master  航迹初始信息
atti_M = sign_set.init_atti_M;    % attitude : roll pitch yaw (deg)
atti_rate_M	= zeros(3,1);   % attitude rate : roll pitch yaw (deg/s)
posi_M		= sign_set.init_posi_M;   % position in GEOREF : lon (deg) lat (deg) alt (m)
veloB_M		= sign_set.init_velo_M;   % velocity on body coordinate : right front up (m/s)
acceB_M		= zeros(3,1);    % acceleration on body coordinate : right front up (m/s^2)
%Rover  航迹初始信息
atti_R = sign_set.init_atti_R;    % attitude : roll pitch yaw (deg)
atti_rate_R	= zeros(3,1);   % attitude rate : roll pitch yaw (deg/s)
posi_R		= sign_set.init_posi_R;   % position in GEOREF : lon (deg) lat (deg) alt (m)
veloB_R		= sign_set.init_velo_R;   % velocity on body coordinate : right front up (m/s)
acceB_R		= zeros(3,1);    % acceleration on body coordinate : right front up (m/s^2)

Xb_real   =zeros(3,sign_set.datalength);

Trace_M=zeros(3,sign_set.datalength);
Trace_R=zeros(3,sign_set.datalength);
Com_Nsat=[];

%% %%%%      读取星历   %%% %%%%
orbit = Ephemeris(sign_set);      % GPS星历

%% %%%%%主函数采样计时开始%%%%%  %%%%%%%
for t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time
     %%    根据航迹求主站、流动站的位置、速度     %%%%
     i=round(1+5*t);
     k=i+1;
     [Xr_M,Vr_M]=user_positions(t,T,atti_M,atti_rate_M,veloB_M,acceB_M,posi_M);
     [Xr_R,Vr_R]=user_positions(t,T,atti_R,atti_rate_R,veloB_R,acceB_R,posi_R);
     Xb_real(:,i)=Xr_M-Xr_R;
     Baseline_real(1,i) = sqrt(Xb_real(1,i)^2+Xb_real(2,i)^2+Xb_real(3,i)^2);          %基线矢量的长度
     Trace_M(:,i)=Xr_M(:,1);                         %  主站航迹
     Trace_R(:,i)=Xr_R(:,1);                            %  流动站航迹
     
     %%       判断主、测站共同可见星   %%%% 
     [com_nsat]=Common_Sat(t, orbit, Xr_M, Xr_R, Deg0);
     Cnsat(:,1)=com_nsat(:,1);
     Cnsat(:,k)=com_nsat(:,2); 
     fprintf('Common Visiable Satellites Number ... ');
     nsat(1,i)=size(Cnsat,1);                % 可见星数目
end   

     %% 求可视卫星到主、测站的传输时间
     [time_M, time_R, ex_M,ey_M,ez_M,ex_R, ey_R, ez_R] = TimeDelay_SR(orbit, Cnsat, Xr_M, Xr_R);
  
     %%  生成频率L1、L2伪距、载波相位观测量
     fprintf('Generating GNSS Observations ... ');
     [pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M,...
         pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R] = GNSS_Observations(time_M, time_R, nsat);
     
     %%  计算载波相位双差观测量     
     [ddpr1, ddpr2, ddph1, ddph2, ddN1, ddN2, Bk3] = GNSS_Carrier_DD(nsat, pr1_M, ph1_M, pr2_M, ph2_M, ...
         N1_M, N2_M, ex_M, ey_M, ez_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R,ex_R, ey_R, ez_R);

%% %%%%%%% Post Processing Mode Calculate Baselinelength%%%%%%%%%%%
% 最小二乘法解基线值，并求模糊度浮点值     
[Xb, ddN1]=LSbaseline_ErrorModel(ddph1,ddN1,ddph2,ddN2, Bk3,Xb0);

ddN0 = ddN1(:,1)*ones(1,size(ddN1,2));  % 基线初始模糊度矩阵
ddN1_dop = ddN1 - ddN0; 
[afixed,t1,Ps1]=LamCacAmb(ddph1,ddN1,Bk3,ddN1_dop,lamada1,Xb);

% 利用Lambda法求解的初始模糊度求取所有模糊度
ddNF1 =afixed(:,1)*ones(1,size(ddN1,2))+ddN1_dop;    % 所有历元模糊度  ???
% 模糊度固定后再求基线向量
[Xb,Q_Xb]=LSbaseline_Ambfixed(ddph1,ddNF1,Bk3,Xb);

for i=1: sign_set.datalength
      Baselinelength(i,1) = sqrt(Xb(1,i)^2+Xb(2,i)^2+Xb(3,i)^2);          %基线矢量的长度
      Baseline_real(i,1)=sqrt(Xb_real(1,i)^2+Xb_real(2,i)^2+Xb_real(3,i)^2)^2;
      Baseline_error(i,1)=Baseline_real(i,1)-Baselinelength(i,1);
end
%}    
     
%% %%%%%%%%%%%%开始绘制图形%%%%%%%%%%%%%%%
fprintf('Now begin to plot figures ... ');
fig_num=0;
%fig_num=fig_num+1;
%figure (fig_num)
t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time;
%plot(t,Baseline_error(:,1),'b');

fig_num=fig_num+1;
figure (fig_num); 
plot3(Trace_M(1,:),Trace_M(2,:),Trace_M(3,:));
fig_num=fig_num+1;
figure (fig_num); 
plot3(Trace_R(1,:),Trace_R(2,:),Trace_R(3,:));
fig_num=fig_num+1;
figure (fig_num); 
plot(t,nsat(1,:));
 