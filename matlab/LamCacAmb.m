 function [afixed,i,Ps] = LamCacAmb(ddph1,ddN1,Bk3,ddN1_dop,bo1,Xb)
%  LAMBDA3法求解整周模糊度
%  input：ddph1  双差载波观测值（周）   （L1 carrier）
%              ddN1    双差模糊度         （L1 carrier）
%              dde1      方向余弦（3维）  （L1 carrier）
%              G_bo1     载波波长
%              Xb     基线长
%  output：afixed    模糊度整数解
%                    i           历元数
%==========================================================================


%% Initialize ======================================================
global sign_set;
%bo1=lamada1;
i = 5;                             % 观测历元数   2历元需同时观测可见卫星数至少4颗
flag = 1;                        % 标志位 表明模糊度多历元解算过程是否结束
ks=length(sign_set.PRNmat); % ks为卫星个数    

Af = Bk3;          %载波的矩阵方程基线改正数系数A 
%while(flag)
%%  构造载波观测方程 L=A*x+B*y  利用LAMBDA法求解模糊度
L   = [ ];
A   = [ ];
B   = [ ];
N1_dop   = [ ];
%De = 2*(ones(ks-1,ks-1) + eye(ks-1));
De=4*((sign_set.Phase*bo1)^2)*eye(ks-1) ;    %  求权矩阵   ？？
Qe = inv(De);                      %求协方差矩阵
P  = [ ];   
% k 根据需要选取多个时刻的数据进行模糊度求解
 for k = 1:i
        L = [L;bo1*(ddph1(:,2*k)-Af(:,:,2*k)*Xb(:,2*k))];   % 存储k个历元载波观测值   5*k  ???
        B = [B;eye(ks-1)];   % 存储k个历元模糊度系数矩阵
        N1_dop = [N1_dop;ddN1_dop(:,2*k)];    % 存储k个历元模糊度值    ???
    
        P1 = [ P, zeros((ks-1)*(k-1),ks-1) ];       %  对应基线协方差矩阵
        P2 = [ zeros(ks-1,(ks-1)*(k-1)), Qe ];    % 对应模糊度协方差矩阵
        P = [ P1; P2];       % 存储k个历元观测量协方差矩阵
 end    
    B = B*bo1;    
    L = L - N1_dop*bo1;
    af1 = inv(B'*P*B)*B'*P*L;
    Q1 = inv(B'*P*B);                %  Q1为模糊度的方差阵
    %ncands = 4;                             % 模糊度候选值个数
%% LAMBDA法求解模糊度
   %[afixed,sqnorms] = lambda2(af1,Qx1,ncands);
    [afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(af1,Q1,4);
%{
    n=0; Q=0;
%% 判断LAMBDA法求解模糊度是否正确，若不正确增加历元  实际接收机如何判断呢
   if(isequal(afixed(:,1),ddN1(:,1)))    % 为什么只比较一个候选值？？
       flag = 0;
       for m=1:ks-1
         if(isequal(afixed(m,1),afixed(m,2))==0)
          Q=Q+abs((af1(m)-afixed(m,1))/(af1(m)-afixed(m,2)));
          n=n+1;
         end   
       end
      rr=1-1/n*Q;
    else
        i = i+ 1;
        i
   end
%}
%% 利用ratio值检验法判断模糊度是否解算正确
%      P1=afixed(:,1)-afx1;
%      P2=afixed(:,2)-afx1;
%      ratio=(P1'*Qx1*P1)/(P2'*Qx1*P2)
%      if(ratio>1.5)
%            flag = 0;
%            afixed(:,1)
%     else
%         i = i+ 1;
%         i
%     end

%% 利用比较方法判断模糊度解算是否正确
% n=0;  % 残差二次型最小和次小值相等的个数
% Q=0;  % 检验值计算的中间值
% for m=1:ks-1
%     if(isequal(afixed(m,1),afixed(m,2))==0)    
%         Q=Q+abs((afx1(m)-afixed(m,1))/(afx1(m)-afixed(m,2)));
%         n=n+1;
%     end
%         
% end
% %计算检验值
% rr=1-1/n*Q
% 
% dd=0;
% if abs(rr)>0.55
%     flag=0;
%     if(rr>0)
%     dd=1;
%     else
%      dd=0;
%     end       
% else
%    i=i+1;
%    
% end

%end

