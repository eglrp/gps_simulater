function [com_nsat]=Common_Sat(t, orbit, Xr_M, Xr_R, Deg0)
%	    �жϿɼ�����Ŀ     �ο�л�ֽ̲�P48
%	Input: 
%         t                              �źŽ���ʱ�䣨�۲�ʱ�̣�
%		  orbit                       ���ǹ������
%         Xr_M                       ��վECEF����
%         Xr_R                        ����վECEF����
%         Deg0                       ��ֹ�߶Ƚ�
%	Output:
%         com_nsat                �۲�ʱ������վ��ͬ����������Ŀ���߶Ƚ�
%      Design by WuLing    2017-12-27
%==========================================================================
global sign_set;
delta_X=zeros(3,1);
nsat0=[];
nsat1=[];
i=round(1+5*t);
% ����������λ�á��ٶ�  
t_orbit = sign_set.t_sate+t;  
q=1;
p=1;
for j=1:length(sign_set.PRNmat) 
         if sum(orbit(j,:))~=0
              [satp,satv]=sate_posivelo(t_orbit,orbit(j,:));
               %�ж�������վ�ɼ�����
              XM_Geo=ecef2geo(Xr_M);            %����ת��ECEF(X,Y, Z)����>GEO(��, ��, ��)      
              XR_Geo=ecef2geo(Xr_R);
              S0=[-sin(XM_Geo(1))  cos(XM_Geo(1))    0;...
                    -sin(XM_Geo(2))*cos(XM_Geo(1))   -sin(XM_Geo(2))*sin(XM_Geo(1))   cos(XM_Geo(2));...
                     cos(XM_Geo(2))*cos(XM_Geo(1))    cos(XM_Geo(2))*sin(XM_Geo(1))   sin(XM_Geo(2)) ];    %����ת������ ECEF����>վ��ϵ
              S1=[-sin(XR_Geo(1))  cos(XR_Geo(1))    0;...
                    -sin(XR_Geo(2))*cos(XR_Geo(1))   -sin(XR_Geo(2))*sin(XR_Geo(1))   cos(XR_Geo(2));...
                     cos(XR_Geo(2))*cos(XR_Geo(1))    cos(XR_Geo(2))*sin(XR_Geo(1))   sin(XR_Geo(2)) ];               
              delta_X0=satp-Xr_M;
              delta_X1=satp-Xr_R;
              Xsp0=S0*delta_X0;                           %������վ������ϵ�µ�����
              Xsp1=S1*delta_X1;  
              % �߶ȽǼ���
              E_deno0(j,i)=Xsp0(1)^2+Xsp0(2)^2;
              E_deno0(j,i)=sqrt(E_deno0(j,i));
              E_deno1(j,i)=Xsp1(1)^2+Xsp1(2)^2;
              E_deno1(j,i)=sqrt(E_deno1(j,i));
                  if   E_deno0(j,i)==0
                       E_rad0(j,i)=pi/2;
                       E_deg0(j,i)=90;
                  else
                       E_rad0(j,i)=atan(Xsp0(3)/E_deno0(j,i));                 %jΪ���Ǳ�ţ�iΪ���ݵ���
                       E_deg0(j,i)=E_rad0(j,i)*180/pi;
                  end 
                  if   E_deno1(j,i)==0
                       E_rad1(j,i)=pi/2;
                       E_deg1(j,i)=90;
                  else
                       E_rad1(j,i)=atan(Xsp1(3)/E_deno1(j,i));                 %jΪ���Ǳ�ţ�iΪ���ݵ���
                       E_deg1(j,i)=E_rad1(j,i)*180/pi;
                  end
         else
              fprintf('error:\n�ļ��в������� %d �����ǵĹ������\n',sign_set.PRNmat(j));
              E_deno0(j,i)=0;
              E_rad0(j,i)=0;
              E_deg0(j,i)=0;
              E_deno1(j,i)=0;
              E_rad1(j,i)=0;
              E_deg1(j,i)=0;
         end
%%%% �߶ȽǱȽϣ��жϿɼ���%%%%%%
      if E_deg0(j,i)>=Deg0
          nsat0(q,1)=j;                           %�洢�ɼ����Ǳ��
          nsat0(q,2)=E_deg0(j,i);             %�洢�ɼ����Ǹ߶Ƚ�
          q=q+1;
      end    
 % �߶ȽǱȽϣ��жϿɼ���
      if E_deg1(j,i)>=Deg0
          nsat1(p,1)=j;                           %�洢�ɼ����Ǳ��
          nsat1(p,2)=E_deg1(j,i);             %�洢�ɼ����Ǹ߶Ƚ�
          p=p+1;
      end
      %%%% �Ƚ���վ�Ͳ�վ�Ĺ������� %%%%%
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
%%%%ѡȡ�߶Ƚ���������Ϊ�ο��ǲ���������һ��%%%%
clear max;   clear row;  
[max0,row]=max(com_sat(:,2)) ;             %��ȡ��߶ȽǼ����Ǳ��
com_sat1(1,:)=com_sat(row,:);                %��������һ��
com_satt=com_sat;
com_satt(row,:)=[];                                  %�ݴ����ɾ���ƶ���
com_sat1=[com_sat1;com_satt];
com_nsat=com_sat1;