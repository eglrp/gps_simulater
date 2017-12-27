function [Xb,Q_xb]=LSbaseline_Ambfixed(ddph1,ddN1,Bk3,Xb)
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

%%  ���̶����ģ���Ȼش����̣����������������
for i=1: sign_set.datalength  
% ��С���˷�������ʸ��
    ks=length(sign_set.PRNmat);
   
    D1 = 2*(ones(ks-1,ks-1) + eye(ks-1))*(sign_set.Phase*sign_set.G_bo1)^2;        %��Ȩ����  GPS-������ P161ҳ  ��ʽ��6.104��
    P1 = inv(D1);      %��Э������� ��С�����е�Ȩ��ֵ

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
     Xb(:,i)=Xb(:,i)+delta_xb;                   %  ��������
     
end