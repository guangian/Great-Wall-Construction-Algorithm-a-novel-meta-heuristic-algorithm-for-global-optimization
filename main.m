% ======================================================================= %
% Great Wall Construction Algorithm: a novel meta-heuristic algorithm for engineer problems
% ----------------------------------------------------------------------- %
%  Programer: Ziyu Guan
%     E-mail: 1909275821@qq.com
% ----------------------------------------------------------------------- %
% Supervisor: Changjiang Ren
%     E-mail: 971932670@qq.com
% ----------------------------------------------------------------------- %
% Main papers:
%             Unpublished
% ======================================================================= %
clc;clear;close;
EngFunction=["Speed Reducer","Tension/compression spring design",...
        "Rolling Element Bearing","Multiple disk clutch brake design problems",...
        "Pressure vessel design","Three-bar truss design problem",...
        "Design of gear train","Cantilever beam","Minimize I-beam vertical deflection",...
        "Tubular column design","Piston lever","Corrugated bulkhead design",...
        "Car side impact design","Design of welded beam design","Reinforced concrete beam design",...
        "Constraint Problem 1","Constraint Problem 2","Constraint Problem 3"]';
% for i=1:length(EngFunction)
for i=1
    disp('-------------------------------------------------------')
    %% Eng
    Fobj=@EngineerProblem;
    Name=EngFunction(i);
    [ub,lb,dim,fobj] = Fobj(Name);%得到上下边界以及维度
    N=100;
    NFE=1000*N;
%     %% CEC-2017
%     Name=['F',num2str(i)];
%     Fobj=@CEC2017;
%     [ub,lb,dim,fobj] = Fobj(Name);%得到上下边界以及维度
%     dim=10;
%     N=50;
%     NFE=10000*dim;
%     %% Unimodal
%     Name=['F',num2str(i)];
%     Fobj=@unimodalVariableDim;
%     [ub,lb,dim,fobj] = Fobj(Name);%得到上下边界以及维度
%     N=50;
%     NFE=1000*N;
    %% 重复次数
    RT=25;
    %% GWCA
    for j=1:RT
        GWCA_Value.N=N;
        GWCA_Value.MaxFe=NFE;
        GWCA_Value.P=7.93;
        GWCA_Value.Q=0.646;
        GWCA_Value.T=2;
        GWCA_Value.m=2.93;
        [Score(j,:),Position,Convergence,PositionValue]=GWCA(GWCA_Value,lb,ub,dim,fobj);
    end
    disp(strcat("GWCA  ",Name,": Best:",string(min(Score)),"; Mean:",string(mean(Score))));
end
