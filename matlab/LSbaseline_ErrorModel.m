function [Xb, ddN1]=LSbaseline_ErrorModel(nsat,ddpr1, ddpr2, ddph1, ddph2, ddN1, ddN2, Bk3)
%	�ز��۲ⷽ��������С���˷����M-R����ʸ��
%	Input: 
%         ddph1 ddph2         L1��L2Ƶ��˫���ز��۲�ֵ
%		  ddN1  ddN2           L1��L2Ƶ��˫��ģ����
%         dde1                       L1  3ά�������Ҿ���(nsat-1)*3  *datalength  
%         Xb0                         α�൥�㶨λ���߳�ֵ
%	Output:
%         Xb                       M-R���߳���
%        ddN1                    L1 ˫��ģ����

%==========================================================================
%% Initialize ======================================================
global sign_set;

%%  ��С����������ʸ��
for i=1: sign_set.datalength  
% ��С���˷�������ʸ��
    ks=nsat(1,i);                % �۲�ʱ�̿ɼ�����Ŀ
   
    %% �������ģ��
    D1 = 2*(ones(ks-1,ks-1) + eye(ks-1))*(sign_set.Phase*sign_set.G_bo1)^2;        %��Ȩ����  GPS-������ P161ҳ  ��ʽ��6.104��
    P1 = inv(D1);      %��Э������� ��С�����е�Ȩ��ֵ

    %% ����˫���ģ�ͣ��۲�ģ�ͣ�
     Af = Bk3(:,:,i);          %�ز��ľ��󷽳̻��߸�����ϵ��A      (A B)*[delta_Xb  delta_ddN]'=L
     Bf=  sign_set.G_bo1*eye(ks-1);            % ˫��ģ���ȸ�����ϵ��B
     Lf = sign_set.G_bo1*(ddph1(:,i)+ddN1(:,i)+Af*Xb(:,i));

     N11=Af'*P1*Af; 
     N12=Af'*P1*Bf;
     N21=Bf'*P1*Af;
     N22=Bf'*P1*Bf;
     U1  =Af'*P1*Lf;
     U2  =Bf'*P1*Lf;
     N_0=N21*inv(N11)*N12;
     P_Q=N22-N_0;
     Segma_A=inv(P_Q);
     delta_ddN=Segma_A*(U2-N21*inv(N11)*U1);      %ģ���ȸ�����
     Q_ddn=Segma_A;     %  ģ����Э�������
     delta_xb=inv(N11)*(U1-N12*delta_ddN);    %���߸�����
     Q_xb=inv(N11)+inv(N11)*N12*Segma_A*N21*inv(N11);     % ����Э�������
     Q_xn=-inv(N11)*N12*Segma_A;     % ����-ģ����Э�������
     Q_nx=-Segma_A*N21*inv(N11);
     Xb(:,i)=Xb0+delta_xb;                   %  ��������
     ddN1(:,i)=ddN1(:,i)+delta_ddN;           %  ģ���ȸ����
   
     %{  
     %% ����LAMBDA�㷨�̶�ģ����
    [afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(delta_ddN,Q_ddn,4);     %��? ����Ԫ�̶�
    
    % �����������������
    delta_xbf=delta_xb(:,i)-Q_xn*Q_ddn*(delta_ddN-afixed);
    Q_xf=Q_xb-Q_xn*P_Q*Q_nx;         %Э����
    Xb(:,i)=Xb0+delta_xbf;
    ddN1(:,i)=ddN1(:,i)+afixed;
   %} 
end