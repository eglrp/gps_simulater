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
sign_set.wie   = 7.2921151467e-5;       %������ת���ٶ�
%sign_set.PI_orbit = 3.1415926535898;         % pi value used for orbit computation

%sign_set.ell_a_GPS = 6378137;                          % GPS (WGS-84)      Ellipsoid semi-major axis [m]
%sign_set.ell_f_GPS = 1/298.257223563;                  % GPS (WGS-84)      Ellipsoid flattening
%sign_set.ell_e_GPS  = sqrt(1-(1-sign_set.ell_f_GPS)^2);    % GPS (WGS-84)      Eccentricity
%sign_set.GM_GPS = 3.986005e14;                        % GPS     Gravitational constant * (mass of Earth) [m^3/s^2]

%%  CONSTELLATION SPECIFIC ================================================
%GPS Ƶ��  [Hz]
sign_set.GPS_L1              =1575.42e6;
sign_set.GPS_L2              =1227.60e6;

% GPS�ز�����  [m]
sign_set.GPS_bo1              =0.1903;
sign_set.GPS_bo2              =0.2442;

%% ȱʡ��������================================
sign_set.GPS_Code           =30;        % �趨������������(m)
sign_set.GPS_Phase         =0.03;      % �趨�ز���λ��������(m)

%%  RINEX�ļ�·��==================================
sign_set.ephfile        =...
    'G:\matlab\GNSS��Ե���_WL\RINEX\brdc1980.17n';        

% ȫ��GPS���Ǳ��
sign_set.PRNmat          =[01 02 03 04 05 06 07 08 09 10 ...
                                         11 12 13 14 15 16 17 18 19 20 ...
                                         21 22 23 24 25 26 27 28 29 30 31 32];    

% ���濪ʼʱ��(��)
sign_set.start_time      = 0;
% �������ʱ��(��) 
sign_set.end_time        =1; 
%%%%%%%%%%%%%%%%%%%%%%%%%
% ��������(������ʱ������
sign_set.TDperiod             = 0.2;
%  �������������
sign_set.datalength=1+(sign_set.end_time - sign_set.start_time)/sign_set.TDperiod;

% �ز���ʼ��λ
sign_set.carrInitPhase   = 0;

%% ��ʼ��������=============================================
% ��̬
sign_set.init_atti_M  =[0 0 0]';
sign_set.init_atti_R  =[0 0 0]';

%��ʼ λ��-GEO  
sign_set.init_posi_M  =[118 32 500]';    %��������ϵ
sign_set.init_posi_R  =[118 32  600]';    %��������ϵ         ���߳���2km??

% �ٶ�
sign_set.init_velo_M  =[0 0 0]';
sign_set.init_velo_R  =[0 10 0]';


% ����λ�á��ٶȼ������ʼʱ��(GPS���룩
sign_set.t_sate=86400;        %�������ļ�

% ���ǹ۲��ֹ�߶Ƚ�(��)
sign_set.cutoff=15             
