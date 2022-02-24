%һ�����㷨
clear
clc
%����
a1=[0.850000000000000;1.14000000000000;1.41000000000000;1.24000000000000;0.980000000000000;0.820000000000000;0.790000000000000;1.06000000000000;1.08000000000000;0.990000000000000];
a2=[0.220000000000000;0.340000000000000;0.770000000000000;0.580000000000000;0.370000000000000;0.390000000000000;0.260000000000000;0.430000000000000;0.820000000000000;0.380000000000000];
Dmax=[195;1040;998;799;448;200;147;687;588;540];
b1=[0.00110000000000000;0.000700000000000000;0.00110000000000000;0.000800000000000000;0.000800000000000000;0.000700000000000000;0.00110000000000000;0.000900000000000000;0.000100000000000000;0.00130000000000000];
b2=[0.000800000000000000;0.000800000000000000;0.000700000000000000;0.00100000000000000;0.000900000000000000;0.000800000000000000;0.00130000000000000;0.00100000000000000;0.00110000000000000;0.000800000000000000];
Smax=[1488;795;652;582;412;670;996;544;210;950];
%��ʼ״̬(�޹���)
S=sdpvar(10,1);
D=sdpvar(10,1);
obj=sum(a1.*D-0.5*b1.*D.^2-a2.*S-0.5*b2.*S.^2);
C1=[0<=S<=Smax,0<=D<=Dmax];
C2=[D==S];%�޹���ģʽ
C=[C1,C2];
ops=sdpsettings('solver','gurobi','gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
solvesdp(C,-obj,ops);
lagrant=dual(C2);%��ʼ�������ճ���(������ֵ)
%����ģʽ(һ�����㷨)
S=double(S);
D=double(D);
%% �������쵼���γ�����ͨ·��һ�����㷨(������1Ϊ�쵼��)
% for t=1:100
%     Link1=0.5*(1-eye(10))/9+0.5*eye(10);Link2=[1;zeros(9,1)];
%     D1=max(min((a1-lagrant)./b1,Dmax),0);
%     S1=max(min((lagrant-a2)./b2,Smax),0);
%     E=D1-S1;
%     lagrant=Link1*lagrant+0.0001*Link2*sum(E);%һ���Ա����ĸ���
% end
%% �������쵼������ȫ���ӵ�һ�����㷨(������1Ϊ�쵼��)
% topology=[1,2;2,1;1,3;3,2;3,4;4,3;3,6;6,3;4,5;5,4;6,7;7,6;5,8;8,5;8,6;8,7;8,9;9,8;8,10;10,8];Link=zeros(10);Link2=[1;zeros(9,1)];
% for i=1:20
%     Link(topology(i,1),topology(i,2))=1;
% end
% Link1=zeros(10);%һ����ϵ��
% for i=1:10
%     for j=1:10
%         if i==j
%             Link1(i,j)=0.5;
%         end
%         if Link(i,j)==1;%i��j�ĸ��ڵ�
%             Link1(i,j)=0.5/sum(Link(i,:));
%         end  
%     end
% end
% for t=1:1000
%     D1=max(min((a1-lagrant)./b1,Dmax),0);
%     S1=max(min((lagrant-a2)./b2,Smax),0);
%     E=D1-S1;
%     lagrant=Link1*lagrant+0.0001*Link2*sum(E);%һ���Ա����ĸ���
% end
%% �������쵼��һ�����㷨(���������е�ͼģ�ͼ���)
topology=[1,2;2,1;1,3;3,2;3,4;4,3;3,7;7,3;4,5;5,6;6,5;7,8;8,7;9,7;6,9;9,6;9,10;10,9];Link=zeros(10);
for i=1:18
    Link(topology(i,1),topology(i,2))=1;
end
w=zeros(10);%һ����ϵ��
v=zeros(10);%������ϵ��
for i=1:10
    for j=1:10
        if i==j
            w(i,j)=0.5;
            v(i,j)=0.5;
        end
        if Link(j,i)==1;%j��i�ĸ��ڵ�
            w(i,j)=0.5/sum(Link(:,i));
        end  
        if Link(i,j)==1%i��j�ĸ��ڵ�
            v(i,j)=0.5/sum(Link(i,:));
        end
    end
end
E_last=zeros(10,1);%��һ�ִεĹ���ƻ�
xigma=zeros(10,1);%������
lagrant_data=[lagrant];
E_data=[];
tic
for t=1:1000
    D1=max(min((a1-lagrant)./b1,Dmax),0);
    S1=max(min((lagrant-a2)./b2,Smax),0);
    E=D1-S1;
    lagrant=w*lagrant+0.0001*xigma;%һ���Ա����ĸ���
    xigma=v'*xigma+(E-E_last);%���·�����
    E_last=E;
    lagrant_data=[lagrant_data,lagrant];
    E_data=[E_data,E];
end
toc
lagrant_data=lagrant_data(:,1:1000);
plot(lagrant_data')
UD=a1.*D1-0.5*b1.*D1.^2-lagrant.*D1;
US=lagrant.*S1-a2.*S1-0.5*b2.*S1.^2;
U=a1.*D1-0.5*b1.*D1.^2-a2.*S1-0.5*b2.*S1.^2;
PD=a1-b1.*D1;
PS=a2+b2.*S1;
% %����ģʽ(��ḣ���������)
% S=sdpvar(10,1);
% D=sdpvar(10,1);
% obj=sum(a1.*D-0.5*b1.*D.^2-a2.*S-0.5*b2.*S.^2);
% C1=[0<=S<=Smax,0<=D<=Dmax];
% C2=[sum(D)==sum(S)];%�޹���ģʽ
% C=[C1,C2];
% ops=sdpsettings('solver','gurobi','gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
% solvesdp(C,-obj,ops);
% lagrant=dual(C2);%������