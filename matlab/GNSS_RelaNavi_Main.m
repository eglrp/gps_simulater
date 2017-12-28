 %	Main funciton of GNSS_RelaNavi_Main
%	Input:
%		void
%	Output:
%		void
%
%	Created by Wu Ling  from NUAA NRC at 2017/10/06
%==========================================================================
%% %%������б������ر����д򿪵�ͼ%%%%%%%%
clear all;
close all;
clc;

%% %%%%%%%%%Get initial settings %%%%%%%%%%%%%%%%%
global sign_set
%  ��ó�ʼ������
sign_set = Inital_Settings();     
lamada1=sign_set.GPS_bo1;
lamada2=sign_set.GPS_bo2;
Deg0=sign_set.cutoff;
%posi_GEO    =zeros(3,sign_set.datalength);      %��������ϵ��λ��
T = sign_set.TDperiod;% ��������(������ʱ������
% Master  ������ʼ��Ϣ
atti_M = sign_set.init_atti_M;    % attitude : roll pitch yaw (deg)
atti_rate_M	= zeros(3,1);   % attitude rate : roll pitch yaw (deg/s)
posi_M		= sign_set.init_posi_M;   % position in GEOREF : lon (deg) lat (deg) alt (m)
veloB_M		= sign_set.init_velo_M;   % velocity on body coordinate : right front up (m/s)
acceB_M		= zeros(3,1);    % acceleration on body coordinate : right front up (m/s^2)
%Rover  ������ʼ��Ϣ
atti_R = sign_set.init_atti_R;    % attitude : roll pitch yaw (deg)
atti_rate_R	= zeros(3,1);   % attitude rate : roll pitch yaw (deg/s)
posi_R		= sign_set.init_posi_R;   % position in GEOREF : lon (deg) lat (deg) alt (m)
veloB_R		= sign_set.init_velo_R;   % velocity on body coordinate : right front up (m/s)
acceB_R		= zeros(3,1);    % acceleration on body coordinate : right front up (m/s^2)

Xb_real   =zeros(3,sign_set.datalength);

Trace_M=zeros(3,sign_set.datalength);
Trace_R=zeros(3,sign_set.datalength);
Com_Nsat=[];

%% %%%%      ��ȡ����   %%% %%%%
orbit = Ephemeris(sign_set);      % GPS����

%% %%%%%������������ʱ��ʼ%%%%%  %%%%%%%
for t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time
     %%    ���ݺ�������վ������վ��λ�á��ٶ�     %%%%
     i=round(1+5*t);
     k=i+1;
     [Xr_M,Vr_M]=user_positions(t,T,atti_M,atti_rate_M,veloB_M,acceB_M,posi_M);
     [Xr_R,Vr_R]=user_positions(t,T,atti_R,atti_rate_R,veloB_R,acceB_R,posi_R);
     Xb_real(:,i)=Xr_M-Xr_R;
     Baseline_real(1,i) = sqrt(Xb_real(1,i)^2+Xb_real(2,i)^2+Xb_real(3,i)^2);          %����ʸ���ĳ���
     Trace_M(:,i)=Xr_M(:,1);                         %  ��վ����
     Trace_R(:,i)=Xr_R(:,1);                            %  ����վ����
     
     %%       �ж�������վ��ͬ�ɼ���   %%%% 
     [com_nsat]=Common_Sat(t, orbit, Xr_M, Xr_R, Deg0);
     Cnsat(:,1)=com_nsat(:,1);
     Cnsat(:,k)=com_nsat(:,2); 
     fprintf('Common Visiable Satellites Number ... ');
     nsat(1,i)=size(Cnsat,1);                % �ɼ�����Ŀ
end   

     %% ��������ǵ�������վ�Ĵ���ʱ��
     [time_M, time_R, ex_M,ey_M,ez_M,ex_R, ey_R, ez_R] = TimeDelay_SR(orbit, Cnsat, Xr_M, Xr_R);
  
     %%  ����Ƶ��L1��L2α�ࡢ�ز���λ�۲���
     fprintf('Generating GNSS Observations ... ');
     [pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M,...
         pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R] = GNSS_Observations(time_M, time_R, nsat);
     
     %%  �����ز���λ˫��۲���     
     [ddpr1, ddpr2, ddph1, ddph2, ddN1, ddN2, Bk3] = GNSS_Carrier_DD(nsat, pr1_M, ph1_M, pr2_M, ph2_M, ...
         N1_M, N2_M, ex_M, ey_M, ez_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R,ex_R, ey_R, ez_R);

%% %%%%%%% Post Processing Mode Calculate Baselinelength%%%%%%%%%%%
% ��С���˷������ֵ������ģ���ȸ���ֵ     
[Xb, ddN1]=LSbaseline_ErrorModel(ddph1,ddN1,ddph2,ddN2, Bk3,Xb0);

ddN0 = ddN1(:,1)*ones(1,size(ddN1,2));  % ���߳�ʼģ���Ⱦ���
ddN1_dop = ddN1 - ddN0; 
[afixed,t1,Ps1]=LamCacAmb(ddph1,ddN1,Bk3,ddN1_dop,lamada1,Xb);

% ����Lambda�����ĳ�ʼģ������ȡ����ģ����
ddNF1 =afixed(:,1)*ones(1,size(ddN1,2))+ddN1_dop;    % ������Ԫģ����  ???
% ģ���ȹ̶��������������
[Xb,Q_Xb]=LSbaseline_Ambfixed(ddph1,ddNF1,Bk3,Xb);

for i=1: sign_set.datalength
      Baselinelength(i,1) = sqrt(Xb(1,i)^2+Xb(2,i)^2+Xb(3,i)^2);          %����ʸ���ĳ���
      Baseline_real(i,1)=sqrt(Xb_real(1,i)^2+Xb_real(2,i)^2+Xb_real(3,i)^2)^2;
      Baseline_error(i,1)=Baseline_real(i,1)-Baselinelength(i,1);
end
%}    
     
%% %%%%%%%%%%%%��ʼ����ͼ��%%%%%%%%%%%%%%%
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
 