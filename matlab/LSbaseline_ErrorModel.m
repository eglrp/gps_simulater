function [Xb, ddN1]=LSbaseline_ErrorModel(nsat,ddpr1, ddpr2, ddph1, ddph2, ddN1, ddN2, Bk3)
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

%%  最小二乘求解基线矢量
for i=1: sign_set.datalength  
% 最小二乘法求解基线矢量
    ks=nsat(1,i);                % 观测时刻可见星数目
   
    %% 设置随机模型
    D1 = 2*(ones(ks-1,ks-1) + eye(ks-1))*(sign_set.Phase*sign_set.G_bo1)^2;        %求权矩阵  GPS-王惠南 P161页  公式（6.104）
    P1 = inv(D1);      %求协方差矩阵 最小二乘中的权重值

    %% 构建双差函数模型（观测模型）
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
     Xb(:,i)=Xb0+delta_xb;                   %  基线向量
     ddN1(:,i)=ddN1(:,i)+delta_ddN;           %  模糊度浮点解
   
     %{  
     %% 利用LAMBDA算法固定模糊度
    [afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(delta_ddN,Q_ddn,4);     %？? 单历元固定
    
    % 再求基线向量改正数
    delta_xbf=delta_xb(:,i)-Q_xn*Q_ddn*(delta_ddN-afixed);
    Q_xf=Q_xb-Q_xn*P_Q*Q_nx;         %协方差
    Xb(:,i)=Xb0+delta_xbf;
    ddN1(:,i)=ddN1(:,i)+afixed;
   %} 
end