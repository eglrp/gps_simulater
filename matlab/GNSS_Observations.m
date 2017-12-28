 function [pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M,  pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R] = GNSS_Observations(TimeDelay_M, TimeDelay_R, nsat)
%   ����α�ࡢ�ز���λ�������չ۲���
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
        i=round(1+5*t);              % �������ݵ���
        nsat0=nsat(1,i);                % �۲�ʱ�̿ɼ�����Ŀ

        for j=1:  nsat0         
                    %% Masterվ�� L1��L2Ƶ��α�ࡢ�ز���λ�۲�ֵ������ģ���ȡ�������Ƶ�ƹ۲���
                   R1_M(j,i)   = sign_set.c*TimeDelay_M(j,i);                            %���ǵ���վ����ʵ����
                   pr1_M(j,i)  =R1_M(j,i) + sign_set.GPS_Code*randn;         %��վα��۲���    ���ף�   ��=r+��
                   pr2_M(j,i)  =R1_M(j,i) + sign_set.GPS_Code*randn;         %��վα��۲���    ���ף�
                   ph1_M(j,i)  = R1_M(j,i)/sign_set.GPS_bo1 + floor(R1_M(j,i)/sign_set.GPS_bo1) + sign_set.GPS_Phase*randn;    %L1Ƶ�� �ز���λ�۲�ֵ  ���ܣ� �˦�=r+N+��
                   ph2_M(j,i)  = R1_M(j,i)/sign_set.GPS_bo2 + floor(R1_M(j,i)/sign_set.GPS_bo2) + sign_set.GPS_Phase*randn;    %L2Ƶ�� �ز���λ�۲�ֵ   ���ܣ�        
                   N1_M(j,i)  =floor(R1_M(j,i)/sign_set.GPS_bo1);         %L1Ƶ�� ����ģ����
                   N2_M(j,i)  =floor(R1_M(j,i)/sign_set.GPS_bo2);         %L2Ƶ�� ����ģ����


                %% Rover ��L1��L2Ƶ��α�ࡢ�ز���λ�۲�ֵ������ģ���ȡ�������Ƶ�ƹ۲���
                   R1_R(j,i)   = sign_set.c*TimeDelay_R(j,i);                           %���ǵ�����վ�ľ���
                   pr1_R(j,i)  =R1_R(j,i) + sign_set.GPS_Code*randn;           %����վα��۲���   
                   pr2_R(j,i)  =R1_R(j,i) + sign_set.GPS_Code*randn;           %����վα��۲���    
                   ph1_R(j,i) =R1_R(j,i)/sign_set.GPS_bo1 + floor(R1_R(j,i)/sign_set.GPS_bo1) + sign_set.GPS_Phase*randn;    %L1Ƶ�� �ز���λ�۲�ֵ
                   ph2_R(j,i) =R1_R(j,i)/sign_set.GPS_bo2 +floor(R1_R(j,i)/sign_set.GPS_bo2) + sign_set.GPS_Phase*randn;    %L2Ƶ�� �ز���λ�۲�ֵ             
                   N1_R(j,i)  =floor(R1_R(j,i)/sign_set.GPS_bo1);         %L1Ƶ�� ����ģ����
                   N2_R(j,i)  =floor(R1_R(j,i)/sign_set.GPS_bo2);         %L2Ƶ�� ����ģ����

        end
end
