function [time_M, time_R, ex_M,ey_M,ez_M,ex_R, ey_R, ez_R] = TimeDelay_SR(orbit, Cnsat, Xr_M, Xr_R)
%	compute satellite positions and transit time
%	Input:  
%           t                        �źŽ���ʱ��
%		    orbit                 ���ǹ����Ϣ
%           Xr_M, Xr_R       ����λ��
%           Cnsat                �ɼ�����
%	Output:  
%           time_M           ���ǵ�����M�Ĵ���ʱ��
%           time_R            ���ǵ�����R�Ĵ���ʱ��
%          ex_M, ey_M, ez_M        ���ǵ�����M��������
%          ex_R, ey_R,  ez_R           ���ǵ�����R��������
%% Initialize ======================================================
global sign_set;

for t=sign_set.start_time: sign_set.TDperiod: sign_set.end_time
    nsat0=size(Cnsat,1);                 % �ɼ�����Ŀ
    i=round(1+5*t);                        % �������ݵ���
%% Master 
traceM.posi_ECEF = Xr_M;
%% Rover 
traceR.posi_ECEF = Xr_R;
%% �����վM������ʱ
Xr1 = traceM.posi_ECEF;
t_orbit = sign_set.t_sate+t;    
        for j=1:nsat0
            % ��������λ���ٶ�
           % j,
            k=Cnsat(j,1);
            [Xs,Vs]=sate_posivelo(t_orbit,orbit(k,:));      %�źŽ���ʱ������λ��
            %% �������㴫����ʱ
            Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
            Tp_old1=0;
            while((Tp1-Tp_old1)>1e-12)
                % ����У��������λ��
                [Xs,Vs]=sate_posivelo(t_orbit-Tp1,orbit(k,:));            %�źŷ���ʱ������λ��
                %������ת����
                % ��t_orbit-Tpʱ�̵�ECEF����ϵת����t_orbitʱ�̵�ECEF����ϵ     
                C=[ cos(sign_set.wie*Tp1)    sin(sign_set.wie*Tp1)    0;...
                    -sin(sign_set.wie*Tp1)    cos(sign_set.wie*Tp1)    0;...
                    0              0              1];
                Xs=C*Xs;
                Vs=C*Vs;
                % ����������ʱ
                Tp_old1=Tp1;
                Tp1=sqrt((Xs-Xr1)'*(Xs-Xr1))/sign_set.c;
            end
            time_M(j,i)=Tp1;
            % ������ǵ���վ�ĵ�λ��������
            ex_M(j,i)= (Xs(1)-Xr1(1))/sqrt((Xs-Xr1)'*(Xs-Xr1));      %  x����
            ey_M(j,i)= (Xs(2)-Xr1(2))/sqrt((Xs-Xr1)'*(Xs-Xr1));      %  y����
            ez_M(j,i)= (Xs(3)-Xr1(3))/sqrt((Xs-Xr1)'*(Xs-Xr1));      %  z����  
        end
        % �������վR�Ĵ�����ʱ
        Xr2 = traceR.posi_ECEF;
        for j=1:nsat0
            % j,
            k=Cnsat(j,1);
            % ��������λ���ٶ�
            [Xs,Vs]=sate_posivelo(t_orbit,orbit(k,:));
            %% �������㴫����ʱ
            Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
            Tp_old2=0;
            while((Tp2-Tp_old2)>1e-12)
                % ����У��������λ��
                [Xs,Vs]=sate_posivelo(t_orbit-Tp2,orbit(k,:));
                % ��t_orbit-Tpʱ�̵�ECEF����ϵת����t_orbitʱ�̵�ECEF����ϵ
                C=[ cos(sign_set.wie*Tp2)    sin(sign_set.wie*Tp2)    0;...
                    -sin(sign_set.wie*Tp2)    cos(sign_set.wie*Tp2)    0;...
                    0              0              1];
                Xs=C*Xs;
                Vs=C*Vs;
                % ����������ʱ
                Tp_old2=Tp2;
                Tp2=sqrt((Xs-Xr2)'*(Xs-Xr2))/sign_set.c;
            end
            time_R(j,i)=Tp2;
            % ������ǵ�����վ�ĵ�λ��������
            ex_R(j,i)= (Xs(1)-Xr2(1))/sqrt((Xs-Xr2)'*(Xs-Xr2));      %  x����
            ey_R(j,i)= (Xs(2)-Xr2(2))/sqrt((Xs-Xr2)'*(Xs-Xr2));      %  y����
            ez_R(j,i)= (Xs(3)-Xr2(3))/sqrt((Xs-Xr2)'*(Xs-Xr2));      %  z����  
        end
end
