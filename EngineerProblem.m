function [ub,lb,dim,fobj] =EngineerProblem(F)
switch F
    %Speed Reducer
    case 'Speed Reducer'
        fobj = @F1;
        lb=[2.6, 0.7, 17, 7.3, 7.3, 2.9, 5];
        ub=[3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5];
        dim=7;
        %Tension/compression spring design
    case 'Tension/compression spring design'
        fobj = @F2;
        lb=[0.05,0.25,2.00];
        ub= [2,1.3,15.0];
        dim=3;
        %Pressure vessel design
    case 'Pressure vessel design'
        fobj = @F3;
        lb=[0.51,0.51,10,10];
        ub=[99.49,99.49,200,200];
        dim=4;
        %Three-bar truss design problem
    case 'Three-bar truss design problem'
        fobj = @F4;
        lb=0*ones(1,2);
        ub=1*ones(1,2);
        dim=2;
        %Design of gear train
    case 'Design of gear train'
        fobj = @F5;
        lb=12*ones(1,4);
        ub=60*ones(1,4);
        dim=4;
        %Cantilever beam
    case 'Cantilever beam'
        fobj = @F6;
        lb=[0.01 0.01 0.01 0.01 0.01];
        ub=[100 100 100 100 100];
        dim=5;
        %Minimize I-beam vertical deflection
    case 'Minimize I-beam vertical deflection'
        fobj = @F7;
        lb=[10 10 0.9 0.9];
        ub=[80 50 5.0 5.0];
        dim=4;
        %Tubular column design
    case 'Tubular column design'
        fobj = @F8;
        lb=[2 0.2];
        ub=[14 0.8];
        dim=2;
        %Piston lever
    case 'Piston lever'
        fobj = @F9;
        lb=[0.05 0.05 0.05 0.05];
        ub=[500 500 500 120];
        dim=4;
        %Corrugated bulkhead design
    case 'Corrugated bulkhead design'
        fobj = @F10;
        lb=[0 0 0 0];
        ub=[100 100 100 5];
        dim=4;
        % Car side impact design
    case 'Car side impact design'
        fobj = @F11;
        lb=[0.50 0.50 0.50 0.50 0.50 0.50 0.50 0 0 -30 -30];
        ub= [1.50 1.50 1.50 1.50 1.50 1.50 1.50 1 1 +30 +30];
        dim=11;
        %Design of welded beam design
    case 'Design of welded beam design'
        fobj = @F12;
        lb=[0.1 0.1 0.1 0.1];
        ub= [2 10 10 2];
        dim=4;
        %Reinforced concrete beam design
    case 'Reinforced concrete beam design'
        fobj = @F13;
        lb=[0 0 5];
        ub=[1 1 10];
        dim=3;
        %Multiple disk clutch brake design problems
    case 'Multiple disk clutch brake design problems'
        fobj=@F14;
        lb=[60 90 1 0 2];
        ub=[80 110 3 1000 9];
        dim=5;
        %Rolling Element Bearing
    case 'Rolling Element Bearing'%結果取负值
        fobj=@F15;
        dim=10;
        D=160;
        d=90;
        lb=[0.5*(D+d) 0.15*(D-d) 4 0.515 0.515 0.4 0.6 0.3 0.02 0.6];
        ub=[0.6*(D+d) 0.45*(D-d) 50 0.6 0.6 0.5 0.7 0.4 0.1 0.85];
        %Constraint Problem 1
    case 'Constraint Problem 1'
        fobj=@F16;
        lb=[-10 -10 -10 -10 -10 -10 -10];
        ub=[10 10 10 10 10 10 10];
        dim=7;
        %Constraint Problem 2
    case 'Constraint Problem 2'
        fobj=@F17;
        lb=[78 33 27 27 27];
        ub=[102 45 45 45 45];
        dim=5;
        %Constraint Problem 3
    case 'Constraint Problem 3'
        fobj=@F18;
        dim=10;
        lb=ones(1,dim)*0;
        ub=ones(1,dim)*1;
        %Constraint Problem 4
    case 'Constraint Problem 4'
        fobj=@F19;
        dim=3;
        lb=ones(1,dim)*0;
        ub=ones(1,dim)*10;
end
end
function [o,Data]=F1(x)
o = 0;
o(1)=0.7854*x(:,1).*x(:,2).^2.*(3.3333.*x(:,3).^2+14.9334.*x(:,3)-43.0934)-1.508.*x(:,1).*(x(:,6).^2+x(:,7).^2).....
    +7.477.*(x(:,6).^3+x(:,7).^3)+0.7854.*(x(:,4).*x(:,6).^2+x(:,5).*x(:,7).^2);
%  o(2)=65856000/(30*10^6*x(4)*x(3)^3);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        g(1)=27/(x(:,1).*x(:,2).^2.*x(:,3))-1;
        g(2)=397.5/(x(:,1).*x(:,2).^2.*x(:,3).^2)-1;
        g(3)=(x(:,4).^(3)*1.93)/(x(:,2).*x(:,6).^4.*x(:,3))-1;
        g(4)=(x(:,5).^3*1.93)/(x(:,2).*x(:,7).^4.*x(:,3))-1;
        g(5)=(sqrt(16.91.*10^6+(745.*x(:,4)./(x(:,2).*x(:,3))).^2))/(110*x(:,6).^(3))-1;
        g(6)=(sqrt(157.5.*10^6+(745.*x(:,5)./(x(:,2).*x(:,3))).^2))/(85*x(:,7).^(3))-1;
        g(7)=(x(:,2).*x(:,3))/(40)-1;
        g(8)=(x(:,2)*5)/(x(:,1))-1;
        g(9)=(x(:,1))/(x(:,2)*12)-1;
        g(10)=(1.5.*x(:,6)+1.9)/(x(:,4))-1;
        g(11)=(1.1.*x(:,7)+1.9)/(x(:,5))-1;
        % 本题无等式约束，故为空；
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % 测试不等式是否成立
        % 不等式的索引函数 H(g)
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % 等式的索引函数
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
%TCSD
function [o,Data]=F2(x)
o = [0];
o(1)=(x(3) + 2)*(x(1)^2)*x(2);
%  o(2)=65856000/(30*10^6*x(4)*x(3)^3);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        g(1)=1-(x(2)^3*x(3))/(71785*(x(1)^4));
        g(2)=((4*x(2)^2 - (x(1)*x(2)))/(12566*(x(2)*x(1)^3 - x(1)^4))) + (1/(5108*x(1)^2))-1;
        g(3)=1 - ((140.45*x(1))/(x(2)^2*x(3)));
        g(4)=((x(2) + x(1))/1.5) - 1;
        % [1-((x(2)^3*x(3))/(71785*x(1)^4));
        %     (4*x(2)^2-x(1)*x(2))/(12566*(x(2)*x(1)^3-x(1)^4))+(1/(5108*x(1)^2))-1;
        %     1-((140.45*x(1))/(x(2)^2*x(3)));
        %     ((x(1)+x(2))/1.5)-1];
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% 压力容器设计问题
%______________________________________________________________________________________
function [o,Data]=F3(x)
o = [0];
o(1)=0.6224*x(1)*x(3)*x(4)+1.7781*x(2)*x(3)^2+3.1661*x(1)^2*x(4)+19.84*x(1)^2*x(3);
%  o(2)=65856000/(30*10^6*x(4)*x(3)^3);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        g(1)=-x(1)+0.0193*x(3);
        g(2)=-x(2)+0.00954*x(3);
        g(3)=-pi*x(3)^2*x(4)-(4/3)*pi*x(3)^3+1296000;
        g(4)=x(4)-240;
        % 本题无等式约束，故为空；
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % 测试不等式是否成立
        % 不等式的索引函数 H(g)
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % 等式的索引函数
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% 三杆桁架设计问题
function [o,Data] = F4(x)
y1 = (((sqrt(2)*x(1)+ x(2))/ ((sqrt(2)*x(1)*x(1))+2*x(1)*x(2)))*2)-2;
y2 = ((( x(2))/ ((sqrt(2)*x(1)*x(1))+2*x(1)*x(2)))*2)-2;
y3= ((1/(sqrt(2)*x(2) + x(1)))*2)-2;
o =(2*sqrt(2)*x(1)+x(2))*100;

if y1 > 0 || y2 > 0 || y3 > 0
    o=o+(abs(y1)+abs(y2)+abs(y3))*1000;
end
Data=x;
end
%Design of gear train
function [o,Data]=F5(x)
x=round(x);
term1=1/6.931;
term2=(x(3)*x(2))/(x(1)*x(4));
o = (term1-term2)^2;
Data=x;
end
% Cantilever beam
function [o,Data]=F6(x)
o = [0];
o(1)=0.0624*(x(1) + x(2) + x(3) + x(4) + x(5));
%  o(2)=65856000/(30*10^6*x(4)*x(3)^3);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        g(1)=(61/(x(1)^3)) + (37/(x(2)^3)) + (19/(x(3)^3)) + (7/(x(4)^3)) + (1/(x(5)^3)) - 1;
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% Minimize I-beam vertical deflection
function [o,Data]=F7(x)
o = [0];
term1 = x(3)*(x(1)-2*x(4))^3/12;
term2 = x(2)*x(4)^3/6;
term3 = 2*x(2)*x(4)*((x(1)-x(4))/2)^2;
o = 5000/(term1+term2+term3);
%  o(2)=65856000/(30*10^6*x(4)*x(3)^3);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        g(1) = 2*x(2)*x(4)+x(3)*(x(1)-2*x(4))-300;
        term1 = x(3)*(x(1)-2*x(4))^3;
        term2 = 2*x(2)*x(4)*(4*x(4)^2+3*x(1)*(x(1)-2*x(4)));
        term3 = (x(1)-2*x(4))*x(3)^3;
        term4 = 2*x(4)*x(2)^3;
        g(2) = ((18*x(1)*10^4)/(term1+term2))+((15*x(2)*10^3)/(term3+term4))-56;
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% Tubular column design
function [o,Data]=F8(x)
o = [0];
o(1)=9.8*x(1)*x(2)+2*x(1);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        g(1)=1.59-x(1)*x(2);
        g(2)=47.4-x(1)*x(2)*(x(1)^2+x(2)^2);
        g(3)=2/x(1)-1;
        g(4)=x(1)/14-1;
        g(5)=0.2/x(2)-1;
        g(6)=x(2)/8-1;
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% Piston lever
function [o,Data]=F9(x)
o = [0];
teta = 0.25*pi; H=x(1); B=x(2); D=x(3); X=x(4);
l2=((X*sin(teta)+H)^2+(B-X*cos(teta))^2)^0.5;
l1=((X-B)^2+H^2)^0.5;
o(1)=0.25*pi*D^2*(l2-l1);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        teta = 0.25*pi; H=x(1); B=x(2); D=x(3); X=x(4); P=1500; Q=10000; L=240; Mmax=1.8e+6;
        R=abs(-X*(X*sin(teta)+H)+H*(B-X*cos(teta)))/sqrt((X-B)^2+H^2);
        F=0.25*pi*P*D^2;
        l2=((X*sin(teta)+H)^2+(B-X*cos(teta))^2)^0.5;
        l1=((X-B)^2+H^2)^0.5;
        g(1)=Q*L*cos(teta)-R*F;
        g(2)=Q*(L-X)-Mmax;
        g(3)=1.2*(l2-l1)-l1;
        g(4)=0.5*D-B;
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% Corrugated bulkhead design
function [o,Data]=F10(x)
o = [0];
b=x(1); h=x(2); l=x(3); t=x(4); ABD = abs(l^2-h^2);
o(1)= (5.885*t*(b+l))/(b+(abs(l^2-h^2))^0.5);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        b=x(1); h=x(2); l=x(3); t=x(4); ABD = abs(l^2-h^2);
        g(1)=-t*h*(0.4*b+l/6)+8.94*(b+(ABD)^0.5);
        g(2)=-t*h^2*(0.2*b+l/12)+2.2*(8.94*(b+(ABD)^0.5))^(4/3);
        g(3)=-t+0.0156*b+0.15;
        g(4)=-t+0.0156*l+0.15;
        g(5)=-t+1.05;
        g(6)=-l+h;
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% Car side impact design
function [o,Data]=F11(x)
o = [0];
% Sections
Sec8 = [0.192 0.345];
Sec9 = [0.192 0.345];
nSec8 = numel(Sec8);
nSec9 = numel(Sec9);
if floor(x(8)*nSec8+1)<=0
    a=1;
else
    a=floor(x(8)*nSec8+1);
end
if floor(x(9)*nSec9+1)<=0
    b=1;
else
    b=floor(x(9)*nSec9+1);
end
x(8) = Sec8(min(a,nSec8));
x(9) = Sec8(min(b,nSec9));
% Objective
o(1)=1.98+4.90*x(1)+6.67*x(2)+6.98*x(3)+4.01*x(4)+1.78*x(5)+2.73*x(7);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^10;
        % Subjective
        Fa =1.16-0.3717*x(2)*x(4)-0.00931*x(2)*x(10)-0.484*x(3)*x(9)+0.01343*x(6)*x(10);
        VCu =0.261-0.0159*x(1)*x(2)-0.188*x(1)*x(8)-0.019*x(2)*x(7)+0.0144*x(3)*x(5)+0.0008757*x(5)*x(10)+0.08045*x(6)*x(9)+0.00139*x(8)*x(11)+0.00001575*x(10)*x(11);
        VCm =0.214+0.00817*x(5)-0.131*x(1)*x(8)-0.0704*x(1)*x(9)+0.03099*x(2)*x(6)-0.018*x(2)*x(7)+0.0208*x(3)*x(8)+0.121*x(3)*x(9)-0.00364*x(5)*x(6)+0.0007715*x(5)*x(10)-0.0005354*x(6)*x(10)+0.00121*x(8)*x(11)+0.00184*x(9)*x(10)-0.02*x(2)^2;
        VCl=0.74-0.61*x(2)-0.163*x(3)*x(8)+0.001232*x(3)*x(10)-0.166*x(7)*x(9)+0.227*x(2)^(2);
        Dur=28.98+3.818*x(3)-4.2*x(1)*x(2)+0.0207*x(5)*x(10)+6.63*x(6)*x(9)-7.7*x(7)*x(8)+0.32*x(9)*x(10);
        Dmr=33.86+2.95*x(3)+0.1792*x(10)-5.057*x(1)*x(2)-11*x(2)*x(8)-0.0215*x(5)*x(10)-9.98*x(7)*x(8)+22*x(8)*x(9);
        Dlr=46.36-9.9*x(2)-12.9*x(1)*x(8)+0.1107*x(3)*x(10);
        Fp=4.72-0.5*x(4)-0.19*x(2)*x(3)-0.0122*x(4)*x(10)+0.009325*x(6)*x(10)+0.000191*x(11)^(2);
        VMBP=10.58-0.674*x(1)*x(2)-1.95*x(2)*x(8)+0.02054*x(3)*x(10)-0.0198*x(4)*x(10)+0.028*x(6)*x(10);
        VFD=16.45-0.489*x(3)*x(7)-0.843*x(5)*x(6)+0.0432*x(9)*x(10)-0.0556*x(9)*x(11)-0.000786*x(11)^(2);
        g = [Fa-1, VCu-0.32, VCm-0.32, VCl-0.32, Dur-32, Dmr-32, Dlr-32, Fp-4, VMBP-9.9, VFD-15.7];
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% Design of welded beam design
function [o,Data]=F12(x)
o = [0];
%
o(1)=1.10471*x(1)^2*x(2)+0.04811*x(3)*x(4)*(14.0+x(2));
%  o(2)=65856000/(30*10^6*x(4)*x(3)^3);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^15;
        Q=6000*(14+x(2)/2);
        D=sqrt(x(2)^2/4+(x(1)+x(3))^2/4);
        J=2*(x(1)*x(2)*sqrt(2)*(x(2)^2/12+(x(1)+x(3))^2/4));
        alpha=6000/(sqrt(2)*x(1)*x(2));
        beta=Q*D/J;
        tau=sqrt(alpha^2+2*alpha*beta*x(2)/(2*D)+beta^2);
        sigma=504000/(x(4)*x(3)^2);
        tmpf=4.013*(30*10^6)/196;
        P=tmpf*sqrt(x(3)^2*x(4)^6/36)*(1-x(3)*sqrt(30/48)/28);
        g(1)=tau-13600;
        g(2)=sigma-30000;
        g(3)=x(1)-x(4);
        g(4)=6000-P;
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end
% Reinforced concrete beam design
function [o,Data]=F13(x)
o = [0];
% Beam Design Params
x1 = [6 6.16 6.32 6.6 7 7.11 7.2 7.8 7.9 8 8.4];
nx1 = numel(x1);
x2 = 28:40;
nx2 = numel(x2);
KK=floor(x(1)*nx1+1);
if KK<=0
    KK=1;
end
As = x1(min(KK,nx1));
KK=floor(x(2)*nx2+1);
if KK<=0
    KK=1;
end
b = x2(min(KK,nx2));
h = x(3);
% Objective
o(1)=29.4*As+0.6*b*h;
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^15;
        KK=floor(x(1)*nx1+1);
        if KK<=0
            KK=1;
        end
        As = x1(min(KK,nx1));
        KK=floor(x(2)*nx2+1);
        if KK<=0
            KK=1;
        end
        b = x2(min(KK,nx2));
        h = x(3);
        g(1)=b/h-4;
        g(2)=180+7.375*As^2/h-As*b;
        % No equality constraint in this problem, so empty;
        geq=[];
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % Test if inequalities hold
        % Index function H(g) for inequalities
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % Index function for equalities
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end

function [o,Data]=F14(x)
% ri=60:80;
% r0=90:110;
% t=[1.5 2 2.5 3];
% f=600:10:1000;
% z=2:9;
x=round(x);
r0=0.0000078;
o=pi*(x(2)^2-x(1)^2)*x(3)*(x(5)+1)*r0;
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        lam=10^10;
        deltar=20;
        lmax=30;
        delta=0.5;
        pmax=1;
        prz=x(4)/(pi*(x(2)^2-x(1)^2));
        n=250;
        vsrmax=10;
        vsr=(2*pi*n*(x(2)^3-x(1)^3))/(90*(x(2)^2-x(1)^2))*0.001;
        iz=55;
        mu=0.5;
        s=1.5;
        ms=40;
        mf=3;
        Tmax=15;
        mh=(2/3)*mu*x(4)*x(5)*((x(2)^3-x(1)^3)/(x(2)^2-x(1)^2))*0.001;
        T=(iz*pi*n)/(30*(mh+mf));
        %-------------------------------
        g=[deltar+x(1)-x(2);
            -lmax+(x(5)+1)*(x(3)+delta);
            prz-pmax;
            prz*vsr-pmax*vsrmax;
            vsr-vsrmax;
            T-Tmax;
            s*ms-mh;
            -T];
        % 本题无等式约束，故为空；
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % 测试不等式是否成立
        % 不等式的索引函数 H(g)
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % 等式的索引函数
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end

function [o,Data]=F15(x)
gama=x(2)/x(1);
fc=37.91*((1+(1.04*((1-gama/1+gama)^1.72)*((x(4)*(2*x(5)-1)/x(5)*(2*x(4)-1))^0.41))^(10/3))^-0.3)*((gama^0.3*(1-gama)^1.39)/(1+gama)^(1/3))*(2*x(4)/(2*x(4)-1))^0.41;

if x(2)<=25.4
    o=-fc*x(3)^(2/3)*x(2)^1.8;
else
    o=-3.647*fc*x(3)^(2/3)*x(2)^1.4;
end
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;lam=10^10;
        D=160;
        d=90;
        bw=30;
        T=D-d-2*x(2);
        phi=2*pi-2*acos(((((D-d)/2)-3*(T/4))^2+(D/2-T/4-x(2))^2-(d/2+T/4)^2)/(2*((D-d)/2-3*(T/4))*(D/2-T/4-x(2))));
        g=[-phi/(2*asin(x(2)/x(1)))+x(3)-1;
            -2*x(2)+x(6)*(D-d);
            -x(7)*(D-d)+2*x(2);
            x(10)*bw-x(2);
            -x(1)+0.5*(D+d);
            -(0.5+x(9))*(D+d)+x(1);
            -0.5*(D-x(1)-x(2))+x(8)*x(2);
            0.515-x(4);
            0.515-x(5)
            x(4)*x(2)-11.033-eps;
            x(5)*x(2)-11.033-eps];
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
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
            if geq==0
                H=0;
            else
                H=1;
            end
        end
    end
Data=x;
end

function [o,Data]=F16(x)
o=(x(1)-10)^2+5*((x(2)-12)^2)+x(3)^4+3*((x(4)-11)^2)+10*(x(5)^6)+7*(x(6)^2)+x(7)^4-4*x(6)*x(7)-10*x(6)-8*x(7);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        %Penalty constant
        lam=10^10;
        g=[-127+2*x(1)^2+3*x(2)^4+x(3)+4*x(4)^2+5*x(5);
            -282+7*x(1)+3*x(2)+10*x(3)^2+x(4)-x(5);
            -196+23*x(1)+x(2)^2+6*x(6)^2-8*x(7);
            4*x(1)^2+x(2)^2-3*x(1)*x(2)+2*x(3)^2+5*x(6)-11*x(7)];
        % 本题无等式约束，故为空；
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % 测试不等式是否成立
        % 不等式的索引函数 H(g)
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % 等式的索引函数
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end

function [o,Data]=F17(x)
o=5.3578547*x(3)^2+0.8356891*x(1)*x(5)+37.293239*x(1)-40792.141;
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        %Penalty constant
        lam=10^10;
        g=[85.334407+0.0056858*x(2)*x(5)+0.0006262*x(1)*x(4)-0.0022053*x(3)*x(5)-92;
            -85.334407-0.0056858*x(2)*x(5)-0.0006262*x(1)*x(4)+0.0022053*x(3)*x(5);
            80.51249+0.0071317*x(2)*x(5)+0.0029955*x(1)*x(2)+0.0021813*x(3)^2-110;
            -80.51249-0.0071317*x(2)*x(5)-0.0029955*x(1)*x(2)-0.0021813*x(3)^2+90;
            9.300961+0.0047026*x(3)*x(5)+0.0012547*x(1)*x(3)+0.0019085*x(3)*x(4)-25;
            -9.300961-0.0047026*x(3)*x(5)-0.0012547*x(1)*x(3)-0.0019085*x(3)*x(4)+20];
        % 本题无等式约束，故为空；
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % 测试不等式是否成立
        % 不等式的索引函数 H(g)
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % 等式的索引函数
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end

function [o,Data]=F18(x)
dim=10;
a=1;
for d=1:dim
    aa=x(d);
    a=aa*a;
end
o=-((sqrt(dim))^dim)*a;
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        %Penalty constant
        lam=10^10;
        nvars=10;
        b=0;
        for i=1:nvars
            bb=x(i)^2;
            b=bb+b;
        end
        g=b-1-eps;
        % 本题无等式约束，故为空；
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % 测试不等式是否成立
        % 不等式的索引函数 H(g)
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % 等式的索引函数
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end

function [o,Data]=F19(x)
o=-((100-(x(1)-5)^2-(x(2)-5)^2-(x(3)-5)^2)/100);
o=o+getnonlinear(x);
    function Z=getnonlinear(x)
        Z=0;
        %Penalty constant
        lam=10^10;
        g=[(x(1)-1)^2+(x(2)-1)^2+(x(3)-1)^2-0.0625;
            (x(1)-2)^2+(x(2)-2)^2+(x(3)-2)^2-0.0625;
            (x(1)-3)^2+(x(2)-3)^2+(x(3)-3)^2-0.0625;
            (x(1)-4)^2+(x(2)-4)^2+(x(3)-4)^2-0.0625;
            (x(1)-5)^2+(x(2)-5)^2+(x(3)-5)^2-0.0625;
            (x(1)-6)^2+(x(2)-6)^2+(x(3)-6)^2-0.0625;
            (x(1)-7)^2+(x(2)-7)^2+(x(3)-7)^2-0.0625;
            (x(1)-8)^2+(x(2)-8)^2+(x(3)-8)^2-0.0625;
            (x(1)-9)^2+(x(2)-9)^2+(x(3)-9)^2-0.0625];
        % 本题无等式约束，故为空；
        geq=[];
        % 应用不等式约束
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % 应用等式约束
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
        % 测试不等式是否成立
        % 不等式的索引函数 H(g)
        function H=getH(g)
            if g<=0
                H=0;
            else
                H=1;
            end
            % 等式的索引函数
            function H=getHeq(geq)
                if geq==0
                    H=0;
                else
                    H=1;
                end
            end
        end
    end
Data=x;
end