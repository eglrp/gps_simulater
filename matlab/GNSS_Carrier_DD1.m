function [ddph1, ddph2, ddN1, ddN2, ddpr1, dde1x,dde1y,dde1z] = GNSS_Carrier_DD1(t,pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_Mx, e1_My,e1_Mz, dop1_M, dop2_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_Rx, e1_Ry, e1_Rz, dop1_R, dop2_R)
%   产生双差载波观测量
%   input：   
%               pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, e1_M, dop1_M, dop2_M   主站观测量        
%               pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, e1_R, dop1_R, dop2_R     流动站观测量  
%   output：    
%               ddph1    ddph2        双差载波观测量(L1 carrier   L2 carrier)   
%               ddN10   ddN20        双差模糊度观测量(L1 carrier   L2 carrier)   
%               ddpr1                   双差伪距观测量
%               dde1                    双差方向余弦矩阵 
%==========================================================================


%% Initialize ======================================================
global sign_set;

nsat=length(sign_set.PRNmat) ;   

%%    产生双差观测量   %%%%%%%%%%%%
for i=1:  nsat         %   i代表哪颗卫星
       for j=1:  sign_set.datalength     %  j代表第几个采样点

           %%单差观测值     单差为站间作差
          dph1(i,j)  =ph1_R(i,j)  - ph1_M(i,j);              %L1频点  载波相位单差观测值
          dph2(i,j)  =ph2_R(i,j)  - ph2_M(i,j);              %L2频点  载波相位单差观测值
          %dph5(i,j)  =ph5_R(i,j)  - ph5_M(i,j);              %L5频点  载波相位单差观测值
           
          dpr1(i,j)=pr1_R(i,j)-pr1_M(i,j);                      % 站间单差伪距观测量
          de1x(i,j)=e1_Mx(i,j)+e1_Rx(i,j);                    % 计算方向余弦 
          de1y(i,j)=e1_My(i,j)+e1_Ry(i,j); 
          de1z(i,j)=e1_Mz(i,j)+e1_Rz(i,j); 
        
          %单差整周模糊度
         dN1(i,j)   = N1_R(i,j)  - N1_M(i,j);              % L1频点 站间单差整周模糊度            
         dN2(i,j)   = N2_R(i,j)  - N2_M(i,j);              % L2频点 站间单差整周模糊度
         %dN5(i,j)   = N5_R(i,j)  - N5_M(i,j);              % L5频点 站间单差整周模糊度      
       end
end

%%  求双差观测值       星间作差
for k=1: nsat-1    %双差为卫星作差, 以第1颗卫星为参考星
       
       %%  双差观测值
       ddph1(k,:)  = dph1(k+1,:)  - dph1(1,:);            %L1双差载波相位
       ddph2(k,:)  = dph2(k+1,:)  - dph2(1,:);            %L2
       %ddph5(k,:)  = dph5(k+1,:)  - dph5(1,:);            %L5

       ddN1(k,:)   = dN1(k+1,:)   - dN1(1,:);  %L1双差模糊度  
       ddN2(k,:)   = dN2(k+1,:)   - dN2(1,:);  %L2双差模糊度
       %ddN5(k,:)   = dN5(k+1,:)   - dN5(1,:);  %L5双差模糊度
       
       ddpr1(k,:) = dpr1(k+1,:)-dpr1(1,:);        % 双差伪距观测量   
       dde1x(k,:) = -0.5*(de1x(k+1,:)-de1x(1,:));        % 双差方向余弦
       dde1y(k,:) = -0.5*(de1y(k+1,:)-de1y(1,:)); 
       dde1z(k,:) = -0.5*(de1z(k+1,:)-de1z(1,:)); 
end
