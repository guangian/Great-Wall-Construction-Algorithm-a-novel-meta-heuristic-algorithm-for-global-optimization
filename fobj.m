function [o,g,h]=fobj(x,func_num)
[f,g,h]=cec20_func(x,func_num);
o=f+getnonlinear(g,h);
end
function Z=getnonlinear(g,geq)
Z=0;
% Penalty constant
lam=10^10;
% 本题不等式约束，故为空；
% 应用不等式约束
for k=1:length(g)
    Z=Z+ lam*g(k)^2*getH(g(k));
end
% 应用等式约束
for k=1:length(geq)
    Z=Z+lam*geq(k)^2*getHeq(geq(k));
end
end
% 测试不等式是否成立
% 不等式的索引函数 H(g)
function H=getH(g)
if g<=0
    H=0;
else
    H=1;
end
end
% 等式的索引函数
function H=getHeq(geq)
if abs(geq)-0.0001<=0
    H=0;
else
    H=1;
end
end