%% Great Wall Construction Algorithm: a novel meta-heuristic algorithm for engineer problems
% Developer and programmer: Ziyu Guan, Changjiang Ren, Peixi Wang
% Contact Info:
% E-Mail:   1909275821@qq.com
% Main paper:
% Great Wall Construction Algorithm
% Journal: Expert Systems With Applications
clc;clear;close all;
warning off
%% 固定随机数种子
noRNG=1;
rng('default')
rng(noRNG)
%% Problem Information
str=["GWCA"];
AlgorithmsInformation=["Algorithms",str];
CEC2017=[1,3:30];
NP=[1:5];
Func=3;%1:Np-hard problems; 2:CEC-2017; 3:Engineer problems
RT=2;%
mm=0;
%% Start
global data
for i=1
    %% Choose problem types
    switch Func
        case 1%Np-hard problems
            FuncName=['F',num2str(NP(i))];
        case 2%CEC-2017
            FuncName=['CEC2017-F',num2str(CEC2017(i))];
        case 3%Engineer problems
            EngFunction=["Speed Reducer","Tension/compression spring design",...
                "Rolling Element Bearing","Multiple disk clutch brake design problems",...
                "Pressure vessel design","Three-bar truss design problem",...
                "Design of gear train","Cantilever beam","Minimize I-beam vertical deflection",...
                "Tubular column design","Piston lever","Corrugated bulkhead design",...
                "Car side impact design","Design of welded beam design","Reinforced concrete beam design",...
                "Constraint Problem 1","Constraint Problem 2","Constraint Problem 3"]';
            GloMin=[2994.4244658,0.012665232788,85538.48,0.313656,6059.714335048436,263.89584338,2.70085714e-12,...
                1.3399576,0.0130741,26.486361473,8.41269832311,6.8429580100808,22.84296954,1.724852308597366,359.2080,...
                680.630057,-30665.5390,-1];
            FuncName=['EngneerProblem-',char(EngFunction(i))];
    end
    [data,option]=GetData(FuncName);
    option.RT=1;%Run times
    option.k=0;%Open Convergence Curve Logging
    option.numAgent=30;        %size of population
    option.maxIteration=1000;    %maximum number of interation
    option.lb=option.lb.*ones(1,option.dim);
    option.ub=option.ub.*ones(1,option.dim);
    ALLPosition=[];
    %% Start calculating GWCA
    for j=1:length(str)
        pause(1e-5)
        Algorithms=@GWCA;%Algorithm
        for k=1:RT
            BestData=Algorithms(option);
            Best_score=BestData.Fitness;
            Best_pos=BestData.BestPosition;
            if option.k==1
                Convergence_curve=BestData.Congervence;
            end
            ALLFitness(k)=Best_score;
            ALLPosition(k,:)=Best_pos;
            if option.k==1
                AllConvergence(:,k)=Convergence_curve;
            end
        end
        if option.k==1
            BestData.AllConvergence=mean(AllConvergence,1);
        end
        BestData.ALLFitness=ALLFitness;
        BestData.StdScore=std(ALLFitness);
        BestData.MeanScore=mean(ALLFitness);
        [BestData.Fitness,index]=min(ALLFitness);
        BestData.BestPosition=ALLPosition(index,:);
        AlgorithmsData{i,j}=BestData;
        OutputData(FuncName,BestData,option,data,char(str(j)))
    end
    FunctionInformation{i}=option;
    FName(i,:)=string(FuncName);
end
if isempty(AlgorithmsData{1})
    AlgorithmsData=AlgorithmsData(i,:);
    FunctionInformation=FunctionInformation(i);
    FName=FName(i);
end
[ConvergencePut,EndOutPut,MeanRankPut,BestSolutionPut]=SortData(AlgorithmsData,str,FName,FunctionInformation,option);
%% Save Data
writematrix(ConvergencePut,'Result\ConvergencePut.xlsx','WriteMode','overwrite')
writematrix(EndOutPut,'Result\Result.xlsx','WriteMode','overwrite')
writematrix(MeanRankPut,'Result\MeanRankPut.xlsx','WriteMode','overwrite')
writematrix(BestSolutionPut,'Result\BestSolutionPut.xlsx','WriteMode','overwrite')