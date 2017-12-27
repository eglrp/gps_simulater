function [Xb,Q_xb]=LSbaseline_Ambfixed(ddph1,ddN1,Bk3,Xb)
%	载波观测方程利用最小二乘法求解M-R基线矢量
%	Input: 
%         ddph1 ddph2         L1、L2频点双差载波观测值
%		  ddN1  ddN2           L1、L2频点双差模糊度
%         dde1                       L1  3维方向余弦矩阵，(nsat-1)*3  *datalength  
%         Xb0                         伪距单点定位基线初值
%	Output:
%         Xb                       M-R基线长度
%        ddN1                    L1 双差模糊度

%==========================================================================
%% Initialize ======================================================
global sign_set;

%%  将固定后的模糊度回代方程，求基线向量改正数
for i=1: sign_set.datalength  
% 最小二乘法求解基线矢量
    ks=length(sign_set.PRNmat);
   
    D1 = 2*(ones(ks-1,ks-1) + eye(ks-1))*(sign_set.Phase*sign_set.G_bo1)^2;        %求权矩阵  GPS-王惠南 P161页  公式（6.104）
    P1 = inv(D1);      %求协方差矩阵 最小二乘中的权重值

     Af = Bk3(:,:,i);          %载波的矩阵方程基线改正数系数A      (A B)*[delta_Xb  delta_ddN]'=L
     Bf=  sign_set.G_bo1*eye(ks-1);            % 双差模糊度改正数系数B
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
     delta_ddN=Segma_A*(U2-N21*inv(N11)*U1);      %模糊度改正数
     Q_ddn=Segma_A;     %  模糊度协方差矩阵
     delta_xb=inv(N11)*(U1-N12*delta_ddN);    %基线改正数
     Q_xb=inv(N11)+inv(N11)*N12*Segma_A*N21*inv(N11);     % 基线协方差矩阵
     Q_xn=-inv(N11)*N12*Segma_A;     % 基线-模糊度协方差矩阵
     Q_nx=-Segma_A*N21*inv(N11);
     Xb(:,i)=Xb(:,i)+delta_xb;                   %  基线向量
     
end