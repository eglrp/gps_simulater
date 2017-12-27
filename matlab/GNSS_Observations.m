 function [pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_M, dop1_M, dop2_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_R, dop1_R, dop2_R, Xb0] = GNSS_Observations(t,XS_tx,VS_tx,Xr_M,Vr_M,Xr_R,Vr_R)
%   产生伪距、载波相位、多普勒观测量
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
for i=1:  nsat         %   i代表哪颗卫星
       for j=1:  sign_set.datalength     %  j代表第几个采样点
            %Master站的 L1、L2、L5频点伪距、载波相位观测值和整周模糊度、多普勒频移观测量
           R1_M(i,j)   = sqrt((XS_tx(i,:)'-Xr_M)'*(XS_tx(i,:)'-Xr_M));                            %卫星到主站的距离
           e1_M(i,:,j)   =(XS_tx(i,:)-Xr_M')/R1_M(i,j);                                               % 卫星到主站的方向余弦，3维矩阵，
           pr1_M(i,j)  =R1_M(i,j) + sign_set.Phase*sign_set.Code*randn;          %主站伪距观测量    ???
           pr2_M(i,j)  =pr1_M(i,j);
           
           ph1_M(i,j)  = R1_M(i,j)/sign_set.G_bo1 - floor(R1_M(i,j)/sign_set.G_bo1) + sign_set.Phase*sign_set.G_bo1*randn;    %L1频点 载波相位观测值
           ph2_M(i,j)  = R1_M(i,j)/sign_set.G_bo2 - floor(R1_M(i,j)/sign_set.G_bo2) + sign_set.Phase*sign_set.G_bo2*randn;    %L2频点 载波相位观测值
           %ph5_M(i,j)  = R1_M(i,j)/sign_set.G_bo5 - floor(R1_M(i,j)/sign_set.G_bo5) + sign_set.Phase*sign_set.G_bo5*randn;    %L5频点 载波相位观测值
                      
           N1_M(i,j)  =floor(R1_M(i,j)/sign_set.G_bo1);         %L1频点 整周模糊度
           N2_M(i,j)  =floor(R1_M(i,j)/sign_set.G_bo2);         %L2频点 整周模糊度
           %N5_M(i,j)  =floor(R1_M(i,j)/sign_set.G_bo5);         %L5频点 整周模糊度          
       
          dop1_M(i,j) = (Vr_M' -VS_tx(i,:))*e1_M(i,:,j)'/sign_set.c*sign_set.L1;       %当卫星与接收机相对远离时，多普勒频移为负，载波相位测量值变大
          dop2_M(i,j) = (Vr_M' -VS_tx(i,:))*e1_M(i,:,j)'/sign_set.c*sign_set.L2;
          %dop5_M(i,j) = (Vr' -VS_tx(i,:))*e1_M(i,j)/sign_set.c*sign_set.L5; 
          e1_Mx(i,j)=(XS_tx(i,1)-Xr_M(1))/R1_M(i,j); 
          e1_My(i,j)=(XS_tx(i,2)-Xr_M(2))/R1_M(i,j);  
          e1_Mz(i,j)=(XS_tx(i,3)-Xr_M(3))/R1_M(i,j); 
          %Rover 的L1、L2、L5频点伪距、载波相位观测值和整周模糊度、多普勒频移观测量
           R1_R(i,j)   = sqrt((XS_tx(i,:)'-Xr_R)'*(XS_tx(i,:)'-Xr_R));                            %卫星到流动站的距离
           e1_R(i,:,j)   =(XS_tx(i,:)-Xr_R')/R1_R(i,j);                                               % 卫星到流动站的方向余弦
           pr1_R(i,j)  =R1_R(i,j) + sign_set.Phase*sign_set.Code*randn;          %流动站伪距观测量    ???
           pr2_R(i,j)  =pr1_R(i,j);
           
           ph1_R(i,j)  = R1_R(i,j)/sign_set.G_bo1 - floor(R1_R(i,j)/sign_set.G_bo1) + sign_set.Phase*sign_set.G_bo1*randn;    %L1频点 载波相位观测值
           ph2_R(i,j)  = R1_R(i,j)/sign_set.G_bo2 - floor(R1_R(i,j)/sign_set.G_bo2) + sign_set.Phase*sign_set.G_bo2*randn;    %L2频点 载波相位观测值
           %ph5_R(i,j)  = R1_R(i,j)/sign_set.G_bo5 - floor(R1_R(i,j)/sign_set.G_bo5) + sign_set.Phase*sign_set.G_bo5*randn;    %L5频点 载波相位观测值
                      
           N1_R(i,j)  =floor(R1_R(i,j)/sign_set.G_bo1);         %L1频点 整周模糊度
           N2_R(i,j)  =floor(R1_R(i,j)/sign_set.G_bo2);         %L2频点 整周模糊度
           %N5_R(i,j)  =floor(R1_R(i,j)/sign_set.G_bo5);         %L5频点 整周模糊度          
       
          dop1_R(i,j) = (Vr_R' -VS_tx(i,:))*e1_R(i,:,j)'/sign_set.c*sign_set.L1;       %当卫星与接收机相对远离时，多普勒频移为负，载波相位测量值变大
          dop2_R(i,j) = (Vr_R' -VS_tx(i,:))*e1_R(i,:,j)'/sign_set.c*sign_set.L2;
          %dop5_R(i,j) = (Vr' -VS_tx(i,:))*e1_R(i,j)/sign_set.c*sign_set.L5; 
          e1_Rx(i,j)=(XS_tx(i,1)-Xr_R(1))/R1_R(i,j); 
          e1_Ry(i,j)=(XS_tx(i,2)-Xr_R(2))/R1_R(i,j);  
          e1_Rz(i,j)=(XS_tx(i,3)-Xr_R(3))/R1_R(i,j); 
       end
end

