%产消者A和产消者B的福利分析，用于分析电能共享市场帕累托改进机理
%% 共享前后电能计划
clear
clc
D1=sdpvar;
S1=sdpvar;
D2=sdpvar;
S2=sdpvar;
obj1=1.1*D1-0.5*0.001*D1^2-0.4*S1-0.5*0.0017*S1^2;
obj2=0.8*D2-0.5*0.0015*D2^2-0.35*S2-0.5*0.0005*S2^2;
obj_total=obj1+obj2;%社会福利
C=[0<=D1<=500,0<=S1<=250,0<=D2<=200,0<=S2<=450];
%C1=[D1==S1,D2==S2];
C1=[D1+D2==S1+S2];
C=[C,C1];
ops=sdpsettings('solver','gurobi');
solvesdp(C,-obj_total,ops);
D1=double(D1);S1=double(S1);D2=double(D2);S2=double(S2);obj1=double(obj1);obj2=double(obj2);obj_total=double(obj_total);
share_price=dual(C1);
surplus1=obj1-share_price*(D1-S1);
surplus2=obj2+share_price*(D1-S1);
%% 共享社会福利
clear
clc
D1_data=[];S1_data=[];D2_data=[];S2_data=[];E_data=[];price_A_data=[];price_B_data=[];surplus1_data=[];surplus2_data=[];
Marginal_D1=[];Marginal_S1=[];Marginal_D2=[];Marginal_S2=[];surplus1_real_data=[];surplus2_real_data=[];
for E=0:0.25:450;
    D1=sdpvar;
    S1=sdpvar;
    D2=sdpvar;
    S2=sdpvar;
    obj1=1.1*D1-0.5*0.001*D1^2-0.4*S1-0.5*0.0017*S1^2;
    obj2=0.8*D2-0.5*0.0015*D2^2-0.35*S2-0.5*0.0005*S2^2;
    C11=[0<=D1<=500,0<=S1<=250];C12=[D1==S1+E];C1=[C11,C12];
    C21=[0<=D2<=200,0<=S2<=450];C22=[D2+E==S2];C2=[C21,C22];
    ops=sdpsettings('solver','gurobi','gurobi.FeasibilityTol',1e-9,'gurobi.OptimalityTol',1e-9);
    solvesdp([C1,C2],-obj1-obj2,ops);
    D1=double(D1);S1=double(S1);D2=double(D2);S2=double(S2);obj1=double(obj1);obj2=double(obj2);
    price_A=dual(C12);price_B=dual(C22);
    surplus1=obj1-price_B*E;surplus2=obj2+price_A*E;%最优反应剩余
    surplus1_real=obj1-0.5*(price_A+price_B)*E;surplus2_real=obj2+0.5*(price_A+price_B)*E;%真实剩余
    D1_data=[D1_data,D1];
    S1_data=[S1_data,S1];
    D2_data=[D2_data,D2];
    S2_data=[S2_data,S2];
    E_data=[E_data,E];
    price_A_data=[price_A_data,price_A];
    price_B_data=[price_B_data,price_B];
    surplus1_data=[surplus1_data,surplus1];
    surplus2_data=[surplus2_data,surplus2];
    surplus1_real_data=[surplus1_real_data,surplus1_real];
    surplus2_real_data=[surplus2_real_data,surplus2_real];
    Marginal_D1=[Marginal_D1,1.1-0.001*D1];
    Marginal_S1=[Marginal_S1,0.4+0.0017*S1];
    Marginal_D2=[Marginal_D2,0.8-0.0015*D2];
    Marginal_S2=[Marginal_S2,0.35+0.0005*S2];
end
welfare_data=1.1*D1_data-0.5*0.001*D1_data.^2-0.4*S1_data-0.5*0.0017*S1_data.^2+0.8*D2_data-0.5*0.0015*D2_data.^2-0.35*S2_data-0.5*0.0005*S2_data.^2;