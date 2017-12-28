function [ddpr1, ddpr2, ddph1, ddph2, ddN1, ddN2, Bk3] = GNSS_Carrier_DD(nsat, pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, Ex_M, Ey_M, Ez_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R,Ex_R, Ey_R, Ez_R)
%   产生双差载波观测量
%   input：   
%               pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, ex_M, ey_M,ez_M   主站观测量        
%               pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, ex_R, ey_R, ez_R     流动站观测量  
%   output：    
%               ddph1    ddph2        双差载波观测量(L1 carrier   L2 carrier)   
%               ddN1   ddN2       双差模糊度观测量(L1 carrier   L2 carrier)   
%               ddpr1, ddpr2                   双差伪距观测量
%               ddex, ddey, ddez                    双差方向余弦矩阵 
%==========================================================================


%% Initialize ======================================================
global sign_set;

for t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time
    i=round(1+5*t);              % 采样数据点数
    nsat0=nsat(1,i)                % 观测时刻可见星数目

        %%    产生双差观测量   %%%%%%%%%%%%
        for j=1:  nsat0         
              %%单差观测值     单差为站间作差
                  dph1(j,i)  =ph1_R(j,i)  - ph1_M(j,i);              %L1频点  载波相位单差观测值
                  dph2(j,i)  =ph2_R(j,i)  - ph2_M(j,i);              %L2频点  载波相位单差观测值        
                  dpr1(j,i)   =pr1_R(j,i)-pr1_M(j,i);                      % 站间单差伪距观测量
                  dpr2(j,i)   =pr2_R(j,i)-pr2_M(j,i);            
                  dex(j,i)    =Ex_M(j,i)-Ex_R(j,i);                       % 计算方向余弦
                  dey(j,i)    =Ey_M(j,i)-Ey_R(j,i);
                  dez(j,i)    =Ez_M(j,i)-Ez_R(j,i);
                  dN1(j,i)   = N1_R(j,i)  - N1_M(j,i);              % L1频点 站间单差整周模糊度            
                  dN2(j,i)   = N2_R(j,i)  - N2_M(j,i);              % L2频点 站间单差整周模糊度
        end

        %%  求双差观测值       星间作差
        for k=1: nsat0-1    %双差为卫星作差, 以第1颗卫星为参考星

               %%  双差观测值
               ddph1(k,i)  = dph1(k+1,i)  - dph1(1,i);            %L1双差载波相位
               ddph2(k,i)  = dph2(k+1,i)  - dph2(1,i);            %L2
               ddN1(k,i)   = dN1(k+1,i)   - dN1(1,i);              %L1双差模糊度  
               ddN2(k,i)   = dN2(k+1,i)   - dN2(1,i);              %L2双差模糊度     
               ddpr1(k,i) = dpr1(k+1,i)-dpr1(1,i);                  % 双差伪距观测量
               ddpr2(k,i) = dpr2(k+1,i)-dpr2(1,i);                  % 双差伪距观测量
               ddex(k,i)  =  dex(k+1,i) - dex(1,i);                   % 双差方向余弦
               ddey(k,i)  =  dey(k+1,i) - dey(1,i);
               ddez(k,i)  =  dez(k+1,i) - dez(1,i);

               Bk3(k,:,i)  =[ddex(k,i),  ddey(k,i), ddez(k,i)];     % 双差方向余弦矩阵  (3维矩阵) ：(nsat0-1)*3  *datalength  
        end
end
