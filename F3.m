function [fit,result,x0]=F3(x)
x0=x;
global data
%%
[~,S]=sort(x);
path0=[data.noS,S(1:data.numLM0),data.noE];
path1=[];
path2=[];
for i=1:length(path0)-1
    [Dist,Path]=graphshortestpath(data.net,path0(i),path0(i+1));
    path1=[path1;path0(i),path0(i+1),Dist]; 
    path2=[path2,Path];
end
fit=sum(path1(:,3));
if nargout>1
    result.fit=fit;    %总目标
    result.path0=path0; %关键路径 （特殊节点）
    result.path1=path1; %关键路径+距离
    result.path2=path2; %详细路径
end
end