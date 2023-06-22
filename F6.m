function [fit,result,x0]=F6(x)
global data
x0=x;
%% 解码
x1=x(1:data.numD);
x(1:data.numD)=[];
[~,S]=sort(x);
selectedC=S(1:data.numSelected);
if isempty(selectedC)
    selectedC=1;
end
selectedC0=data.noC(selectedC);
%% 安排节点的配送路径
[~,S]=sort(x1);
Load=selectedC*0+data.maxLoad;
flagC=selectedC*0;
recording=[];
demand=data.demand(:,1:data.numP);
demandC=zeros(length(selectedC),3);
for i=1:data.numD
    noD=data.noD(S(i));
    noD0=S(i);  
    for noP=1:data.numP
        if demand(noD0,noP)>0
            position=find(Load>demand(noD0,noP));
            if ~isempty(position)
                [D1,no]=min(data.D1(noD,selectedC(position)));
                Load(position(no))=Load(position(no))-demand(noD0,noP);
                noC=selectedC(position(no));
                noC0=selectedC0(position(no));
                noP0=data.noP(noP);
                D2=data.D2(noC,noP)/1000;
                demandC(position(no),noP)=demandC(position(no),noP)+demand(noD0,noP);
                recording=[recording;noP,noC,noD0,noP0,noC0,noD,demand(noD0,noP),D1,D2];
                % 1生产地独立编号 2物流中心独立编号 3需求地独立编号
                % 4生产地统一编号 5物流中心统一编号 6需求地统一编号
                % 7需求 8距离1-需求地到物流中心 9距离2-物流中心到生产地
                demand(noD0,noP)=0;
            end
        end
    end
end
%% 固定成本
C1=sum(data.node(selectedC0,4));
%% 运输成本
C21=data.ct2*sum(recording(:,7).*(recording(:,8)));
C22=data.ct1*sum(recording(:,7).*(recording(:,9)));
%% 可变成本
C3=sum(sum(data.demand.^data.alpha))*data.cb;
%% 库存成本
C4=sum(sum(data.demand))*data.ck/12;
%% 惩罚项-是否所有的需求均被满足
punishiment=sum(sum(demand(demand>0)));
fit=C1+C21+C22+C3+C4+punishiment*1e6;
if nargout>1
    result.fit=fit; %总目标
    result.recording=recording; %详细记录
    % 1生产地独立编号 2物流中心独立编号 3需求地独立编号
    % 4生产地统一编号 5物流中心统一编号 6需求地统一编号
    % 7需求 8距离1-需求地到物流中心 9距离2-物流中心到生产地
    result.selectedC0=selectedC0; %2物流中心独立编号
    result.selectedC=selectedC;   %5物流中心统一编号
    result.C1=C1;
    result.C21=C21;
    result.C22=C22;
    result.C3=C3;
    result.C4=C4;
    result.demandC=demandC; %每个中转的负载
    result.punishiment=punishiment; %多少需求未被满足
end
end