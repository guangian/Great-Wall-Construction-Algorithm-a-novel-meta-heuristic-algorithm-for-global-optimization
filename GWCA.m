function [DestinationFitness,BestPosition,Conver,PositionValue]=GWCA(Value,LB,UB,dim,FCN)
%重力加速度 g
%石块质量 m
%分形维数 C
%坡度 sitar
%淘汰概率 p
%伽马分布参数P
%伽马分布参数Q
%工具磨损程度SL
%% 基本参数
N=Value.N;
MaxIter=Value.MaxFe/N;
P=Value.P;%伽马分布参数P
Q=Value.Q;%伽马分布参数Q
T=Value.T;%工具产生的推力
m=Value.m;%石块质量
%% 定义算法初始参数
Fitness=zeros(N,1);%对适应度值初始值赋0
record_fitness=zeros(1,MaxIter);
g=9.8;%重力加速度
Cmax=exp(3);%分形维数
Cmin=exp(2);
LNP=ceil(N*0.2);%淘汰人数
%% 初始化种群位置和最佳适应度值
Po=init(N,dim,UB,LB);
BestFitness=inf;
Xbest=zeros(1,dim);%对全局极值进行赋0
for i=1:N
    [Fitness(i,:),Po(i,:)]=FCN(Po(i,:));
    if Fitness(i)<BestFitness
        BestFitness=Fitness(i);
        Xbest=Po(i,:);
    end
end
%% 开始迭代循环
for it=1:MaxIter
    for i=1:N
        Index=randi([1,3]);
        if Index==1
            %% 工程师的速度和位置更新机制
            C=log((Cmax-Cmin)*(MaxIter-it)/MaxIter+Cmin);
            H=(1-it/MaxIter);
            sitar=80*rand(1,dim);
            TL=1-it/MaxIter+eps;
            a=(T.*TL)./m-g.*(H./sin(sitar));
            v=a.*C.*gampdf(it,P,Q);
            Study=(-1)^randi([0,1])*(Xbest-Po(i,:)).*rand(1,dim);
            Po(i,:)=Xbest+Study+Po(i,:).*v;                             %Eq.(3)
        elseif Index==2
            %% 力工的更新速度和位置
            C=log((Cmax-Cmin)*(MaxIter-it)/MaxIter+Cmin);
            H=(1-it/MaxIter)+eps;
            sitar=80*rand(1,dim);
            v=m.*g.*(H./sin(sitar)).*C.*gampdf(it/MaxIter,P,Q);
            Index=1:N;
            Index(Index==i)=[];
            a=Fitness(Index)-Fitness(i);
            [~,ide]=min(abs(a));
            improve=sign(Fitness(ide)-Fitness(i))*(Po(ide,:)-Po(i,:)).*v;
            study=(Xbest-Po(i,:)).*rand(1,dim);
            Po(i,:)=Po(i,:)+improve+study;
        elseif Index==3
            %% 奴隶(军人)的速度和位置更新机制
            Po(i,:)=Po(i,:).*gampdf(it/MaxIter,P,Q);
        end
        Flag4ub=Po(i,:)>UB;
        Flag4lb=Po(i,:)<LB;
        Po(i,:)=(Po(i,:).*(~(Flag4ub+Flag4lb)))+UB.*Flag4ub+LB.*Flag4lb;
        [Fitness(i,:),Po(i,:)]=FCN(Po(i,:));
        if Fitness(i,:)<BestFitness
            BestFitness=Fitness(i,:);
            Xbest=Po(i,:);
        end
    end
    %% 人员淘汰机制
    [~,Index]=sort(Fitness);
%     Po=Po(Index,:);
%     if i>=(N-LNP)
    Po(Index(N-LNP+1:end),:)=(UB-LB).*rand(LNP,dim)+LB;
%     Po(N-LNP+1:end,:)=(UB-LB).*rand(LNP,dim)+LB;
%         [Fitness(i,:),Po(i,:)]=FCN(Po(i,:));
%     end
%     [~,Index]=sort(Fitness);
%     Po=Po(Index,:);
    record_fitness(it)=BestFitness;
    PositionValue{1,it}=Xbest;%储存最佳参数
    PositionValue{2,it}=BestFitness;%储存最佳适应度值
    PositionValue{3,it}=Po;%储存种群参数
    PositionValue{4,it}=Fitness;%储存种群适应度值
end
%% 记录每次得到的适应度值曲线和最优参数还有最佳适应度值
Conver=record_fitness;
DestinationFitness=BestFitness;
BestPosition=Xbest;
end
function Y=gampdf(x,a,b)
Y = 1/(gamma(x)*b^a)*x^(a-1)*exp(-x/b);
end