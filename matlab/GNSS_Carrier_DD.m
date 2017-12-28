function [ddpr1, ddpr2, ddph1, ddph2, ddN1, ddN2, Bk3] = GNSS_Carrier_DD(nsat, pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, Ex_M, Ey_M, Ez_M, pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R,Ex_R, Ey_R, Ez_R)
%   ����˫���ز��۲���
%   input��   
%               pr1_M, ph1_M, pr2_M, ph2_M, N1_M, N2_M, ex_M, ey_M,ez_M   ��վ�۲���        
%               pr1_R, ph1_R, pr2_R, ph2_R, N1_R, N2_R, ex_R, ey_R, ez_R     ����վ�۲���  
%   output��    
%               ddph1    ddph2        ˫���ز��۲���(L1 carrier   L2 carrier)   
%               ddN1   ddN2       ˫��ģ���ȹ۲���(L1 carrier   L2 carrier)   
%               ddpr1, ddpr2                   ˫��α��۲���
%               ddex, ddey, ddez                    ˫������Ҿ��� 
%==========================================================================


%% Initialize ======================================================
global sign_set;

for t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time
    i=round(1+5*t);              % �������ݵ���
    nsat0=nsat(1,i)                % �۲�ʱ�̿ɼ�����Ŀ

        %%    ����˫��۲���   %%%%%%%%%%%%
        for j=1:  nsat0         
              %%����۲�ֵ     ����Ϊվ������
                  dph1(j,i)  =ph1_R(j,i)  - ph1_M(j,i);              %L1Ƶ��  �ز���λ����۲�ֵ
                  dph2(j,i)  =ph2_R(j,i)  - ph2_M(j,i);              %L2Ƶ��  �ز���λ����۲�ֵ        
                  dpr1(j,i)   =pr1_R(j,i)-pr1_M(j,i);                      % վ�䵥��α��۲���
                  dpr2(j,i)   =pr2_R(j,i)-pr2_M(j,i);            
                  dex(j,i)    =Ex_M(j,i)-Ex_R(j,i);                       % ���㷽������
                  dey(j,i)    =Ey_M(j,i)-Ey_R(j,i);
                  dez(j,i)    =Ez_M(j,i)-Ez_R(j,i);
                  dN1(j,i)   = N1_R(j,i)  - N1_M(j,i);              % L1Ƶ�� վ�䵥������ģ����            
                  dN2(j,i)   = N2_R(j,i)  - N2_M(j,i);              % L2Ƶ�� վ�䵥������ģ����
        end

        %%  ��˫��۲�ֵ       �Ǽ�����
        for k=1: nsat0-1    %˫��Ϊ��������, �Ե�1������Ϊ�ο���

               %%  ˫��۲�ֵ
               ddph1(k,i)  = dph1(k+1,i)  - dph1(1,i);            %L1˫���ز���λ
               ddph2(k,i)  = dph2(k+1,i)  - dph2(1,i);            %L2
               ddN1(k,i)   = dN1(k+1,i)   - dN1(1,i);              %L1˫��ģ����  
               ddN2(k,i)   = dN2(k+1,i)   - dN2(1,i);              %L2˫��ģ����     
               ddpr1(k,i) = dpr1(k+1,i)-dpr1(1,i);                  % ˫��α��۲���
               ddpr2(k,i) = dpr2(k+1,i)-dpr2(1,i);                  % ˫��α��۲���
               ddex(k,i)  =  dex(k+1,i) - dex(1,i);                   % ˫�������
               ddey(k,i)  =  dey(k+1,i) - dey(1,i);
               ddez(k,i)  =  dez(k+1,i) - dez(1,i);

               Bk3(k,:,i)  =[ddex(k,i),  ddey(k,i), ddez(k,i)];     % ˫������Ҿ���  (3ά����) ��(nsat0-1)*3  *datalength  
        end
end
