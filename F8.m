function [fit,result,x0]=F8(x)
x0=x;
global data
%%
typeL=round(x.*(data.numLine-1))+1;
R=data.Line(typeL,1).*data.distance(:,2)/data.base;
X=data.Line(typeL,2).*data.distance(:,2)/data.base;
I=data.Line(typeL,3).*data.distance(:,2);
C=data.Line(typeL,4).*data.distance(:,2);
data.mpc.branch(:,3)=R;
data.mpc.branch(:,4)=X;
[MVAbase, bus, gen, gencost, branch, f, success, et]  = ...
                runopf(data.mpc);
lamdaP=bus(:,14);
P=bus(:,3);
[loss, fchg, tchg] = get_losses(MVAbase, bus, branch);
if success==1
    punishiment=0;
else
    punishiment=10000;
end
%%
A=sum(lamdaP.*P);
B=sum(real(loss))*data.lamdaLoss;
C=sum(C);
fit=(A+B)*8760*15+C+punishiment;
%% 
if nargout>1
    result.fit=fit;
    result.bus=bus;
    result.branch=branch;
    result.loss=loss;
    result.lamdaP=lamdaP;
    result.P=P;
    result.A=A; %发电成本
    result.B=B; %损耗
    result.C=C;  %线路成本
end
end