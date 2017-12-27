 function [pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_M, dop1_M, dop2_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_R, dop1_R, dop2_R, Xb0] = GNSS_Observations(t,XS_tx,VS_tx,Xr_M,Vr_M,Xr_R,Vr_R)
%   ����α�ࡢ�ز���λ�������չ۲���
% INPUT:
%   Xs = satellite positions (X,Y,Z)
%   Vs =satellite velocity 
%   Xr_M = Master real position (X,Y,Z)
%   Vr_M = Master real velocity
%   Xr_R = Rover real position (X,Y,Z)
%   Vr_R = Rover real velocity
% OUTPUT: 
%   pr1_M = code observation (L1 carrier)
%   ph1_M = phase observation (L1 carrier)
%   pr2_M = code observation (L2 carrier)
%   ph2_M = phase observation (L2 carrier)
%   N1_M  = ambiguity (L1 carrier)
%   N2_M  = ambiguity (L2 carrier)
%   e1_M   = direction cosine matrix
%   dop1_M = Doppler observation (L1 carrier)
%   dop2_M = Doppler observation (L2 carrier)
%   pr1_R = code observation (L1 carrier)
%   ph1_R = phase observation (L1 carrier)
%   pr2_R = code observation (L2 carrier)
%   ph2_R = phase observation (L2 carrier)
%   N1_R  = ambiguity (L1 carrier)
%   N2_R  = ambiguity (L2 carrier)
%   e1_R   = direction cosine matrix
%   dop1_R = Doppler observation (L1 carrier)
%   dop2_R = Doppler observation (L2 carrier)         
%   Xb0  =  initial baseline vector
%==========================================================================


%% Initialize ======================================================
global sign_set;

Xb0=Xr_M-Xr_R;
nsat=length(sign_set.PRNmat) ;     
for i=1:  nsat         %   i�����Ŀ�����
       for j=1:  sign_set.datalength     %  j����ڼ���������
            %Masterվ�� L1��L2��L5Ƶ��α�ࡢ�ز���λ�۲�ֵ������ģ���ȡ�������Ƶ�ƹ۲���
           R1_M(i,j)   = sqrt((XS_tx(i,:)'-Xr_M)'*(XS_tx(i,:)'-Xr_M));                            %���ǵ���վ�ľ���
           e1_M(i,:,j)   =(XS_tx(i,:)-Xr_M')/R1_M(i,j);                                               % ���ǵ���վ�ķ������ң�3ά����
           pr1_M(i,j)  =R1_M(i,j) + sign_set.Phase*sign_set.Code*randn;          %��վα��۲���    ???
           pr2_M(i,j)  =pr1_M(i,j);
           
           ph1_M(i,j)  = R1_M(i,j)/sign_set.G_bo1 - floor(R1_M(i,j)/sign_set.G_bo1) + sign_set.Phase*sign_set.G_bo1*randn;    %L1Ƶ�� �ز���λ�۲�ֵ
           ph2_M(i,j)  = R1_M(i,j)/sign_set.G_bo2 - floor(R1_M(i,j)/sign_set.G_bo2) + sign_set.Phase*sign_set.G_bo2*randn;    %L2Ƶ�� �ز���λ�۲�ֵ
           %ph5_M(i,j)  = R1_M(i,j)/sign_set.G_bo5 - floor(R1_M(i,j)/sign_set.G_bo5) + sign_set.Phase*sign_set.G_bo5*randn;    %L5Ƶ�� �ز���λ�۲�ֵ
                      
           N1_M(i,j)  =floor(R1_M(i,j)/sign_set.G_bo1);         %L1Ƶ�� ����ģ����
           N2_M(i,j)  =floor(R1_M(i,j)/sign_set.G_bo2);         %L2Ƶ�� ����ģ����
           %N5_M(i,j)  =floor(R1_M(i,j)/sign_set.G_bo5);         %L5Ƶ�� ����ģ����          
       
          dop1_M(i,j) = (Vr_M' -VS_tx(i,:))*e1_M(i,:,j)'/sign_set.c*sign_set.L1;       %����������ջ����Զ��ʱ��������Ƶ��Ϊ�����ز���λ����ֵ���
          dop2_M(i,j) = (Vr_M' -VS_tx(i,:))*e1_M(i,:,j)'/sign_set.c*sign_set.L2;
          %dop5_M(i,j) = (Vr' -VS_tx(i,:))*e1_M(i,j)/sign_set.c*sign_set.L5; 
          e1_Mx(i,j)=(XS_tx(i,1)-Xr_M(1))/R1_M(i,j); 
          e1_My(i,j)=(XS_tx(i,2)-Xr_M(2))/R1_M(i,j);  
          e1_Mz(i,j)=(XS_tx(i,3)-Xr_M(3))/R1_M(i,j); 
          %Rover ��L1��L2��L5Ƶ��α�ࡢ�ز���λ�۲�ֵ������ģ���ȡ�������Ƶ�ƹ۲���
           R1_R(i,j)   = sqrt((XS_tx(i,:)'-Xr_R)'*(XS_tx(i,:)'-Xr_R));                            %���ǵ�����վ�ľ���
           e1_R(i,:,j)   =(XS_tx(i,:)-Xr_R')/R1_R(i,j);                                               % ���ǵ�����վ�ķ�������
           pr1_R(i,j)  =R1_R(i,j) + sign_set.Phase*sign_set.Code*randn;          %����վα��۲���    ???
           pr2_R(i,j)  =pr1_R(i,j);
           
           ph1_R(i,j)  = R1_R(i,j)/sign_set.G_bo1 - floor(R1_R(i,j)/sign_set.G_bo1) + sign_set.Phase*sign_set.G_bo1*randn;    %L1Ƶ�� �ز���λ�۲�ֵ
           ph2_R(i,j)  = R1_R(i,j)/sign_set.G_bo2 - floor(R1_R(i,j)/sign_set.G_bo2) + sign_set.Phase*sign_set.G_bo2*randn;    %L2Ƶ�� �ز���λ�۲�ֵ
           %ph5_R(i,j)  = R1_R(i,j)/sign_set.G_bo5 - floor(R1_R(i,j)/sign_set.G_bo5) + sign_set.Phase*sign_set.G_bo5*randn;    %L5Ƶ�� �ز���λ�۲�ֵ
                      
           N1_R(i,j)  =floor(R1_R(i,j)/sign_set.G_bo1);         %L1Ƶ�� ����ģ����
           N2_R(i,j)  =floor(R1_R(i,j)/sign_set.G_bo2);         %L2Ƶ�� ����ģ����
           %N5_R(i,j)  =floor(R1_R(i,j)/sign_set.G_bo5);         %L5Ƶ�� ����ģ����          
       
          dop1_R(i,j) = (Vr_R' -VS_tx(i,:))*e1_R(i,:,j)'/sign_set.c*sign_set.L1;       %����������ջ����Զ��ʱ��������Ƶ��Ϊ�����ز���λ����ֵ���
          dop2_R(i,j) = (Vr_R' -VS_tx(i,:))*e1_R(i,:,j)'/sign_set.c*sign_set.L2;
          %dop5_R(i,j) = (Vr' -VS_tx(i,:))*e1_R(i,j)/sign_set.c*sign_set.L5; 
          e1_Rx(i,j)=(XS_tx(i,1)-Xr_R(1))/R1_R(i,j); 
          e1_Ry(i,j)=(XS_tx(i,2)-Xr_R(2))/R1_R(i,j);  
          e1_Rz(i,j)=(XS_tx(i,3)-Xr_R(3))/R1_R(i,j); 
       end
end

