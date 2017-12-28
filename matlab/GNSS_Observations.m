 function [pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M,  pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R] = GNSS_Observations(TimeDelay_M, TimeDelay_R, nsat)
%   产生伪距、载波相位、多普勒观测量
% INPUT:
%   time_M = Master time delay
%   time_R = Rover time delay
% OUTPUT: 
%   pr1_M = code observation (L1 carrier)
%   ph1_M = phase observation (L1 carrier)
%   pr2_M = code observation (L2 carrier)
%   ph2_M = phase observation (L2 carrier)
%   N1_M  = ambiguity (L1 carrier)
%   N2_M  = ambiguity (L2 carrier)
%   pr1_R = code observation (L1 carrier)
%   ph1_R = phase observation (L1 carrier)  
%   pr2_R = code observation (L2 carrier)
%   ph2_R = phase observation (L2 carrier)
%   N1_R  = ambiguity (L1 carrier)
%   N2_R  = ambiguity (L2 carrier)

%==========================================================================


%% Initialize ======================================================
global sign_set;

for t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time
        i=round(1+5*t);              % 采样数据点数
        nsat0=nsat(1,i);                % 观测时刻可见星数目

        for j=1:  nsat0         
                    %% Master站的 L1、L2频点伪距、载波相位观测值和整周模糊度、多普勒频移观测量
                   R1_M(j,i)   = sign_set.c*TimeDelay_M(j,i);                            %卫星到主站的真实距离
                   pr1_M(j,i)  =R1_M(j,i) + sign_set.GPS_Code*randn;         %主站伪距观测量    （米）   ρ=r+ε
                   pr2_M(j,i)  =R1_M(j,i) + sign_set.GPS_Code*randn;         %主站伪距观测量    （米）
                   ph1_M(j,i)  = R1_M(j,i)/sign_set.GPS_bo1 + floor(R1_M(j,i)/sign_set.GPS_bo1) + sign_set.GPS_Phase*randn;    %L1频点 载波相位观测值  （周） λφ=r+N+ε
                   ph2_M(j,i)  = R1_M(j,i)/sign_set.GPS_bo2 + floor(R1_M(j,i)/sign_set.GPS_bo2) + sign_set.GPS_Phase*randn;    %L2频点 载波相位观测值   （周）        
                   N1_M(j,i)  =floor(R1_M(j,i)/sign_set.GPS_bo1);         %L1频点 整周模糊度
                   N2_M(j,i)  =floor(R1_M(j,i)/sign_set.GPS_bo2);         %L2频点 整周模糊度


                %% Rover 的L1、L2频点伪距、载波相位观测值和整周模糊度、多普勒频移观测量
                   R1_R(j,i)   = sign_set.c*TimeDelay_R(j,i);                           %卫星到流动站的距离
                   pr1_R(j,i)  =R1_R(j,i) + sign_set.GPS_Code*randn;           %流动站伪距观测量   
                   pr2_R(j,i)  =R1_R(j,i) + sign_set.GPS_Code*randn;           %流动站伪距观测量    
                   ph1_R(j,i) =R1_R(j,i)/sign_set.GPS_bo1 + floor(R1_R(j,i)/sign_set.GPS_bo1) + sign_set.GPS_Phase*randn;    %L1频点 载波相位观测值
                   ph2_R(j,i) =R1_R(j,i)/sign_set.GPS_bo2 +floor(R1_R(j,i)/sign_set.GPS_bo2) + sign_set.GPS_Phase*randn;    %L2频点 载波相位观测值             
                   N1_R(j,i)  =floor(R1_R(j,i)/sign_set.GPS_bo1);         %L1频点 整周模糊度
                   N2_R(j,i)  =floor(R1_R(j,i)/sign_set.GPS_bo2);         %L2频点 整周模糊度

        end
end
