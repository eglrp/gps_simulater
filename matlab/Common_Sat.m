function [com_nsat]=Common_Sat(t, orbit, Xr_M, Xr_R, Deg0)
%	    判断可见星数目     参考谢钢教材P48
%	Input: 
%         t                              信号接收时间（观测时刻）
%		  orbit                       卫星轨道参数
%         Xr_M                       主站ECEF坐标
%         Xr_R                        流动站ECEF坐标
%         Deg0                       截止高度角
%	Output:
%         com_nsat                观测时刻主测站共同可视卫星数目及高度角
%      Design by WuLing    2017-12-27
%==========================================================================
global sign_set;
delta_X=zeros(3,1);
nsat0=[];
nsat1=[];
i=round(1+5*t);
% 求所有卫星位置、速度  
t_orbit = sign_set.t_sate+t;  
q=1;
p=1;
for j=1:length(sign_set.PRNmat) 
         if sum(orbit(j,:))~=0
              [satp,satv]=sate_posivelo(t_orbit,orbit(j,:));
               %判断主、测站可见卫星
              XM_Geo=ecef2geo(Xr_M);            %坐标转换ECEF(X,Y, Z)——>GEO(λ, Φ, Η)      
              XR_Geo=ecef2geo(Xr_R);
              S0=[-sin(XM_Geo(1))  cos(XM_Geo(1))    0;...
                    -sin(XM_Geo(2))*cos(XM_Geo(1))   -sin(XM_Geo(2))*sin(XM_Geo(1))   cos(XM_Geo(2));...
                     cos(XM_Geo(2))*cos(XM_Geo(1))    cos(XM_Geo(2))*sin(XM_Geo(1))   sin(XM_Geo(2)) ];    %坐标转换矩阵 ECEF——>站心系
              S1=[-sin(XR_Geo(1))  cos(XR_Geo(1))    0;...
                    -sin(XR_Geo(2))*cos(XR_Geo(1))   -sin(XR_Geo(2))*sin(XR_Geo(1))   cos(XR_Geo(2));...
                     cos(XR_Geo(2))*cos(XR_Geo(1))    cos(XR_Geo(2))*sin(XR_Geo(1))   sin(XR_Geo(2)) ];               
              delta_X0=satp-Xr_M;
              delta_X1=satp-Xr_R;
              Xsp0=S0*delta_X0;                           %卫星在站心坐标系下的坐标
              Xsp1=S1*delta_X1;  
              % 高度角计算
              E_deno0(j,i)=Xsp0(1)^2+Xsp0(2)^2;
              E_deno0(j,i)=sqrt(E_deno0(j,i));
              E_deno1(j,i)=Xsp1(1)^2+Xsp1(2)^2;
              E_deno1(j,i)=sqrt(E_deno1(j,i));
                  if   E_deno0(j,i)==0
                       E_rad0(j,i)=pi/2;
                       E_deg0(j,i)=90;
                  else
                       E_rad0(j,i)=atan(Xsp0(3)/E_deno0(j,i));                 %j为卫星编号，i为数据点数
                       E_deg0(j,i)=E_rad0(j,i)*180/pi;
                  end 
                  if   E_deno1(j,i)==0
                       E_rad1(j,i)=pi/2;
                       E_deg1(j,i)=90;
                  else
                       E_rad1(j,i)=atan(Xsp1(3)/E_deno1(j,i));                 %j为卫星编号，i为数据点数
                       E_deg1(j,i)=E_rad1(j,i)*180/pi;
                  end
         else
              fprintf('error:\n文件中不包含第 %d 颗卫星的轨道参数\n',sign_set.PRNmat(j));
              E_deno0(j,i)=0;
              E_rad0(j,i)=0;
              E_deg0(j,i)=0;
              E_deno1(j,i)=0;
              E_rad1(j,i)=0;
              E_deg1(j,i)=0;
         end
%%%% 高度角比较，判断可见星%%%%%%
      if E_deg0(j,i)>=Deg0
          nsat0(q,1)=j;                           %存储可见卫星编号
          nsat0(q,2)=E_deg0(j,i);             %存储可见卫星高度角
          q=q+1;
      end    
 % 高度角比较，判断可见星
      if E_deg1(j,i)>=Deg0
          nsat1(p,1)=j;                           %存储可见卫星编号
          nsat1(p,2)=E_deg1(j,i);             %存储可见卫星高度角
          p=p+1;
      end
      %%%% 比较主站和测站的共视卫星 %%%%%
        dimen_com=0;
        for dimen_0=1:q-1
            for dimen=1:p-1
                if nsat0(dimen_0,1)==nsat1(dimen,1)
                    dimen_com=dimen_com+1;
                    com_sat(dimen_com,1)=nsat0(dimen_0,1);
                    com_sat(dimen_com,2)=nsat0(dimen_0,2);
                end
            end 
        end            
end
%%%%选取高度角最大的卫星为参考星并排序至第一行%%%%
clear max;   clear row;  
[max0,row]=max(com_sat(:,2)) ;             %提取最高度角及卫星编号
com_sat1(1,:)=com_sat(row,:);                %排序至第一颗
com_satt=com_sat;
com_satt(row,:)=[];                                  %暂存矩阵删除移动行
com_sat1=[com_sat1;com_satt];
com_nsat=com_sat1;