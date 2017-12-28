function sign_set = Inital_Settings()
%	Initial settings of the GNSS_Rela_Navi_Sim
%	Input:
%		void
%	Output:
%		settings	- simulation settings (structure)
%
%	Created by Wu Ling from NUAA NRC at 2017/09/30
%==========================================================================
%%  Constellation  setting =========================================================
sign_set.c = 299792458;        % light speed (m/s)
sign_set.g = 9.7803698;         % gravity (m/s^2)
sign_set.wie   = 7.2921151467e-5;       %地球自转角速度
%sign_set.PI_orbit = 3.1415926535898;         % pi value used for orbit computation

%sign_set.ell_a_GPS = 6378137;                          % GPS (WGS-84)      Ellipsoid semi-major axis [m]
%sign_set.ell_f_GPS = 1/298.257223563;                  % GPS (WGS-84)      Ellipsoid flattening
%sign_set.ell_e_GPS  = sqrt(1-(1-sign_set.ell_f_GPS)^2);    % GPS (WGS-84)      Eccentricity
%sign_set.GM_GPS = 3.986005e14;                        % GPS     Gravitational constant * (mass of Earth) [m^3/s^2]

%%  CONSTELLATION SPECIFIC ================================================
%GPS 频率  [Hz]
sign_set.GPS_L1              =1575.42e6;
sign_set.GPS_L2              =1227.60e6;

% GPS载波波长  [m]
sign_set.GPS_bo1              =0.1903;
sign_set.GPS_bo2              =0.2442;

%% 缺省精度设置================================
sign_set.GPS_Code           =30;        % 设定测距码测量精度(m)
sign_set.GPS_Phase         =0.03;      % 设定载波相位测量精度(m)

%%  RINEX文件路径==================================
sign_set.ephfile        =...
    'G:\matlab\GNSS相对导航_WL\RINEX\brdc1980.17n';        

% 全部GPS卫星编号
sign_set.PRNmat          =[01 02 03 04 05 06 07 08 09 10 ...
                                         11 12 13 14 15 16 17 18 19 20 ...
                                         21 22 23 24 25 26 27 28 29 30 31 32];    

% 仿真开始时间(秒)
sign_set.start_time      = 0;
% 仿真结束时间(秒) 
sign_set.end_time        =1; 
%%%%%%%%%%%%%%%%%%%%%%%%%
% 航迹生成(传输延时）周期
sign_set.TDperiod             = 0.2;
%  航迹采样点个数
sign_set.datalength=1+(sign_set.end_time - sign_set.start_time)/sign_set.TDperiod;

% 载波初始相位
sign_set.carrInitPhase   = 0;

%% 初始航迹设置=============================================
% 姿态
sign_set.init_atti_M  =[0 0 0]';
sign_set.init_atti_R  =[0 0 0]';

%初始 位置-GEO  
sign_set.init_posi_M  =[118 32 500]';    %地理坐标系
sign_set.init_posi_R  =[118 32  600]';    %地理坐标系         基线长度2km??

% 速度
sign_set.init_velo_M  =[0 0 0]';
sign_set.init_velo_R  =[0 10 0]';


% 卫星位置、速度计算的起始时间(GPS周秒）
sign_set.t_sate=86400;        %由星历文件

% 卫星观测截止高度角(度)
sign_set.cutoff=15             
