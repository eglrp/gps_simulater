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
posi_GEO    =zeros(3,sign_set.datalength);      %��������ϵ��λ��
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
%Master��L1��L2Ƶ����ز���λֵ
ph1_M          =zeros(length(sign_set.PRNmat),sign_set.datalength);
ph2_M          =zeros(length(sign_set.PRNmat),sign_set.datalength);

%Rover��L1��L2Ƶ����ز���λֵ
ph1_R      =zeros(length(sign_set.PRNmat),sign_set.datalength);
ph2_R      =zeros(length(sign_set.PRNmat),sign_set.datalength);

%���ǵ�����˫�������ʸ�� 
Bk3   =zeros(length(sign_set.PRNmat)-1,3,sign_set.datalength);
%˫������ģ����
ddN1  = zeros(length(sign_set.PRNmat)-1,sign_set.datalength);
ddN2  = zeros(length(sign_set.PRNmat)-1,sign_set.datalength);
Xb_real   =zeros(3,sign_set.datalength);

Trace_M=zeros(3,sign_set.datalength);
Trace_R=zeros(3,sign_set.datalength);
Com_Nsat=[];
%% %%%%%%%%%* �����������α�ࡢ�ز���λ�۲���*%%%%%%%%%%%%%
%fprintf('Generating GNSS Observations ... ');
%% %%%%      ��ȡ����   %%% %%%%
orbit = Ephemeris(sign_set);      % GPS����

%% %%%%%������������ʱ��ʼ%%%%%  %%%%%%%
for t=sign_set.start_time: T: sign_set.end_time
     %%%%    ���ݺ�������վ������վ��λ�á��ٶ�     %%%%
     i=round(1+5*t);
     k=i+1;
     [Xr_M,Vr_M]=user_positions(t,T,atti_M,atti_rate_M,veloB_M,acceB_M,posi_M);
     [Xr_R,Vr_R]=user_positions(t,T,atti_R,atti_rate_R,veloB_R,acceB_R,posi_R);
     Xb_real(:,i)=Xr_M-Xr_R;
     Baseline_real(1,i) = sqrt(Xb_real(1,i)^2+Xb_real(2,i)^2+Xb_real(3,i)^2);          %����ʸ���ĳ���
     Trace_M(:,i)=Xr_M(:,1);                         %  ��վ����
     Trace_R(:,i)=Xr_R(:,1);                            %  ����վ����

     %%%%       �ж�������վ��ͬ�ɼ���   %%%% 
     [com_nsat]=Common_Sat(t, orbit, Xr_M, Xr_R, Deg0);
     Cnsat(:,1)=com_nsat(:,1);
     Cnsat(:,k)=com_nsat(:,2); 
     nsat(1,i)=size(Cnsat,1);                 % �ɼ�����Ŀ
     % ��������ǵ�������վ�Ĵ���ʱ��
     [time_M, time_R] = TimeDelay_SR(t, orbit, Cnsat, Xr_M, Xr_R);
     %  ���ɹ۲�������
     
     
     
end



     
     

%{
     Xb_real(:,i)=Xr_M-Xr_R;
     %  ���ɹ۲�������
     [pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_Mx, e1_My, e1_Mz, dop1_M, dop2_M, ...
         pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_Rx, e1_Ry, e1_Rz, dop1_R, dop2_R, Xb0] = GNSS_Observations1(t,XS_tx,VS_tx,Xr_M,Vr_M,Xr_R,Vr_R);
     %  �����ز���λ˫��۲���
      [ddph1, ddph2, ddN1, ddN2, ddpr1, dde1x,dde1y,dde1z] = GNSS_Carrier_DD1(t,...
          pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_Mx, e1_My, e1_Mz, dop1_M, dop2_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_Rx, e1_Ry, e1_Rz, dop1_R, dop2_R);
%end
disp('done.');
for i=1: sign_set.datalength
    for k=1:length(sign_set.PRNmat)-1
        Bk3(k,:,i) =[dde1x(k,i),dde1y(k,i),dde1z(k,i)];    
    end 
end
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
 