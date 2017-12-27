 function [afixed,i,Ps] = LamCacAmb(ddph1,ddN1,Bk3,ddN1_dop,bo1,Xb)
%  LAMBDA3���������ģ����
%  input��ddph1  ˫���ز��۲�ֵ���ܣ�   ��L1 carrier��
%              ddN1    ˫��ģ����         ��L1 carrier��
%              dde1      �������ң�3ά��  ��L1 carrier��
%              G_bo1     �ز�����
%              Xb     ���߳�
%  output��afixed    ģ����������
%                    i           ��Ԫ��
%==========================================================================


%% Initialize ======================================================
global sign_set;
%bo1=lamada1;
i = 5;                             % �۲���Ԫ��   2��Ԫ��ͬʱ�۲�ɼ�����������4��
flag = 1;                        % ��־λ ����ģ���ȶ���Ԫ��������Ƿ����
ks=length(sign_set.PRNmat); % ksΪ���Ǹ���    

Af = Bk3;          %�ز��ľ��󷽳̻��߸�����ϵ��A 
%while(flag)
%%  �����ز��۲ⷽ�� L=A*x+B*y  ����LAMBDA�����ģ����
L   = [ ];
A   = [ ];
B   = [ ];
N1_dop   = [ ];
%De = 2*(ones(ks-1,ks-1) + eye(ks-1));
De=4*((sign_set.Phase*bo1)^2)*eye(ks-1) ;    %  ��Ȩ����   ����
Qe = inv(De);                      %��Э�������
P  = [ ];   
% k ������Ҫѡȡ���ʱ�̵����ݽ���ģ�������
 for k = 1:i
        L = [L;bo1*(ddph1(:,2*k)-Af(:,:,2*k)*Xb(:,2*k))];   % �洢k����Ԫ�ز��۲�ֵ   5*k  ???
        B = [B;eye(ks-1)];   % �洢k����Ԫģ����ϵ������
        N1_dop = [N1_dop;ddN1_dop(:,2*k)];    % �洢k����Ԫģ����ֵ    ???
    
        P1 = [ P, zeros((ks-1)*(k-1),ks-1) ];       %  ��Ӧ����Э�������
        P2 = [ zeros(ks-1,(ks-1)*(k-1)), Qe ];    % ��Ӧģ����Э�������
        P = [ P1; P2];       % �洢k����Ԫ�۲���Э�������
 end    
    B = B*bo1;    
    L = L - N1_dop*bo1;
    af1 = inv(B'*P*B)*B'*P*L;
    Q1 = inv(B'*P*B);                %  Q1Ϊģ���ȵķ�����
    %ncands = 4;                             % ģ���Ⱥ�ѡֵ����
%% LAMBDA�����ģ����
   %[afixed,sqnorms] = lambda2(af1,Qx1,ncands);
    [afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(af1,Q1,4);
%{
    n=0; Q=0;
%% �ж�LAMBDA�����ģ�����Ƿ���ȷ��������ȷ������Ԫ  ʵ�ʽ��ջ�����ж���
   if(isequal(afixed(:,1),ddN1(:,1)))    % Ϊʲôֻ�Ƚ�һ����ѡֵ����
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
%% ����ratioֵ���鷨�ж�ģ�����Ƿ������ȷ
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

%% ���ñȽϷ����ж�ģ���Ƚ����Ƿ���ȷ
% n=0;  % �в��������С�ʹ�Сֵ��ȵĸ���
% Q=0;  % ����ֵ������м�ֵ
% for m=1:ks-1
%     if(isequal(afixed(m,1),afixed(m,2))==0)    
%         Q=Q+abs((afx1(m)-afixed(m,1))/(afx1(m)-afixed(m,2)));
%         n=n+1;
%     end
%         
% end
% %�������ֵ
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

