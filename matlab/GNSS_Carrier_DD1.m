function [ddph1, ddph2, ddN1, ddN2, ddpr1, dde1x,dde1y,dde1z] = GNSS_Carrier_DD1(t,pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_Mx, e1_My,e1_Mz, dop1_M, dop2_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_Rx, e1_Ry, e1_Rz, dop1_R, dop2_R)
%   ����˫���ز��۲���
%   input��   
%               pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_M, dop1_M, dop2_M   ��վ�۲���        
%               pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_R, dop1_R, dop2_R     ����վ�۲���  
%   output��    
%               ddph1    ddph2        ˫���ز��۲���(L1 carrier   L2 carrier)   
%               ddN10   ddN20        ˫��ģ���ȹ۲���(L1 carrier   L2 carrier)   
%               ddpr1                   ˫��α��۲���
%               dde1                    ˫������Ҿ��� 
%==========================================================================


%% Initialize ======================================================
global sign_set;

nsat=length(sign_set.PRNmat) ;   

%%    ����˫��۲���   %%%%%%%%%%%%
for i=1:  nsat         %   i�����Ŀ�����
       for j=1:  sign_set.datalength     %  j����ڼ���������

           %%����۲�ֵ     ����Ϊվ������
          dph1(i,j)  =ph1_R(i,j)  - ph1_M(i,j);              %L1Ƶ��  �ز���λ����۲�ֵ
          dph2(i,j)  =ph2_R(i,j)  - ph2_M(i,j);              %L2Ƶ��  �ز���λ����۲�ֵ
          %dph5(i,j)  =ph5_R(i,j)  - ph5_M(i,j);              %L5Ƶ��  �ز���λ����۲�ֵ
           
          dpr1(i,j)=pr1_R(i,j)-pr1_M(i,j);                      % վ�䵥��α��۲���
          de1x(i,j)=e1_Mx(i,j)+e1_Rx(i,j);                    % ���㷽������ 
          de1y(i,j)=e1_My(i,j)+e1_Ry(i,j); 
          de1z(i,j)=e1_Mz(i,j)+e1_Rz(i,j); 
        
          %��������ģ����
         dN1(i,j)   = N1_R(i,j)  - N1_M(i,j);              % L1Ƶ�� վ�䵥������ģ����            
         dN2(i,j)   = N2_R(i,j)  - N2_M(i,j);              % L2Ƶ�� վ�䵥������ģ����
         %dN5(i,j)   = N5_R(i,j)  - N5_M(i,j);              % L5Ƶ�� վ�䵥������ģ����      
       end
end

%%  ��˫��۲�ֵ       �Ǽ�����
for k=1: nsat-1    %˫��Ϊ��������, �Ե�1������Ϊ�ο���
       
       %%  ˫��۲�ֵ
       ddph1(k,:)  = dph1(k+1,:)  - dph1(1,:);            %L1˫���ز���λ
       ddph2(k,:)  = dph2(k+1,:)  - dph2(1,:);            %L2
       %ddph5(k,:)  = dph5(k+1,:)  - dph5(1,:);            %L5

       ddN1(k,:)   = dN1(k+1,:)   - dN1(1,:);  %L1˫��ģ����  
       ddN2(k,:)   = dN2(k+1,:)   - dN2(1,:);  %L2˫��ģ����
       %ddN5(k,:)   = dN5(k+1,:)   - dN5(1,:);  %L5˫��ģ����
       
       ddpr1(k,:) = dpr1(k+1,:)-dpr1(1,:);        % ˫��α��۲���   
       dde1x(k,:) = -0.5*(de1x(k+1,:)-de1x(1,:));        % ˫�������
       dde1y(k,:) = -0.5*(de1y(k+1,:)-de1y(1,:)); 
       dde1z(k,:) = -0.5*(de1z(k+1,:)-de1z(1,:)); 
end
