function BestData=GWCA(option)
N=option.numAgent;
MaxIter=option.maxIteration;
LB=option.lb;
UB=option.ub;
dim=option.dim;
FCN=option.fobj;
tic;
RT=option.RT;
for irun=1:RT
    % A memory for counting number of function evaluations
    neval = 0;
    Worker1=zeros(1,dim);   Worker1_fit=inf;
    Worker2=zeros(1,dim);   Worker2_fit=inf;
    Worker3=zeros(1,dim);   Worker3_fit=inf;
    %% Initializing
    Fitness=zeros(N,1);
    SL=1;
    T=8.3;
    g=9.8;
    m=3;
    e=0.1;
    P=9;
    Q=6;
    Cmax=exp(3);
    Cmin=exp(2);
    LNP=ceil(N*e);
    %% Initialize the population position and the best fitness value
    Po=initial(N,dim,UB,LB);                           %Eq.(1)
    for i=1:N
        Fitness(i,:)=FCN(Po(i,:));
        if Fitness(i)<Worker1_fit
            Worker1_fit=Fitness(i);  Worker1=Po(i,:);
        elseif Fitness(i)>Worker1_fit && Fitness(i)<Worker2_fit
            Worker2_fit=Fitness(i);  Worker2=Po(i,:);
        elseif Fitness(i)>Worker1_fit && Fitness(i)>Worker2_fit && Fitness(i)<Worker3_fit
            Worker3_fit=Fitness(i);  Worker3=Po(i,:);
        end
    end
    Pbest=Fitness;
    P_Pobest=Po;
    %% Start
    for it=1:MaxIter
        for i=1:N
            Index=randi([1,3]);
            if Index==1
                %% engineer movement
                C=log((Cmax-Cmin)*(MaxIter-it)/MaxIter+Cmin);              %Eq.(7)
                H=(1-it/MaxIter);                                          %Eq.(6)
                sitar=80*rand(1,dim);
                TL=SL-it/MaxIter+eps;
                a=(T.*TL)./m-g.*(H./sin(sitar));
                v=a.*C.*gampdf(it,P,Q);                                    %Eq.(5)
                Study=(-1)^randi([0,1])*(Worker1-Po(i,:)).*rand(1,dim);%Eq.(4)
                Po(i,:)=Worker1+Study+Po(i,:).*v.*rand(1,dim);                            %Eq.(3)
            elseif Index==2
                %% soldier movement
                C=log((Cmax-Cmin)*(MaxIter-it)/MaxIter+Cmin);
                H=(1-it/MaxIter)+eps;
                sitar=80*rand(1,dim);
                v=m.*g.*(H./sin(sitar)).*C.*gampdf(it/MaxIter,P,Q);        %Eq.(9)
                Index=1:N;
                Index(Index==i)=[];
                a=Fitness(Index)-Fitness(i);
                [~,ide]=min(abs(a));
                improve=sign(Fitness(ide)-Fitness(i))*(Po(ide,:)-Po(i,:)).*v.*rand(1,dim);
                study=(Worker2-Po(i,:)).*rand(1,dim);
                Po(i,:)=Po(i,:)+improve+study;                             %Eq.(10)
            elseif Index==3
                %% labor movement
                Po(i,:)=Po(i,:)+2*(Worker3-Po(i,:)).*rand(1,dim)+(P_Pobest(i,:)-Po(i,:)).*gampdf(it/MaxIter,P,Q);                   %Eq.(11)
            end
            %% constraint bound
            Flag4ub=Po(i,:)>UB;
            Flag4lb=Po(i,:)<LB;
            Po(i,:)=(Po(i,:).*(~(Flag4ub+Flag4lb)))+UB.*Flag4ub+LB.*Flag4lb;
            Fitness(i,:)=FCN(Po(i,:));
            neval = neval + 1;
            %update local optimum
            if Fitness(i,:)<Pbest(i,:)
                Pbest(i,:)=Fitness(i,:);
                P_Pobest(i,:)=Po(i,:);
            end
            %update leader
            if Fitness(i)<Worker1_fit
                Worker1_fit=Fitness(i);  Worker1=Po(i,:);
            elseif Fitness(i)>Worker1_fit && Fitness(i)<Worker2_fit
                Worker2_fit=Fitness(i);  Worker2=Po(i,:);
            elseif Fitness(i)>Worker1_fit && Fitness(i)>Worker2_fit && Fitness(i)<Worker3_fit
                Worker3_fit=Fitness(i);  Worker3=Po(i,:);
            end
        end
        %% Personnel elimination mechanism
        [~,Index]=sort(Fitness,'descend');
        Po(Index(1:LNP),:)=(UB-LB).*rand(LNP,dim)+LB;
        %% sort fitness
        record_fitness(it)=Worker1_fit;
        disp(['Iter',num2str(it),'   Best Fitness',num2str(Worker1_fit)])
    end
    %% Record the fitness value curve and optimal parameters obtained each time, as well as the best fitness value
    Nfe(irun,:) = neval;
    Conver(irun,:)=record_fitness;
    DestinationFitness(irun,:)=Worker1_fit;
    BestPosition(irun,:)=Worker1;
end
BestData.Congervence = mean(record_fitness,1);
BestData.MeanNfe = mean(Nfe);
BestData.Conver=mean(Conver,1);
BestData.ALLFitness=DestinationFitness;
BestData.StdScore=std(DestinationFitness);
BestData.MeanScore=mean(DestinationFitness);
[BestData.Fitness,MINindex]=min(DestinationFitness);
BestData.BestPosition=BestPosition(MINindex,:);
end
function Y=gampdf(x,a,b)
Y = 1/(gamma(x)*b^a)*x^(a-1)*exp(-x/b);
end
%% Level 1: Initializing
function Positions=initial(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries
cxl=rand(SearchAgents_no,dim);
for j=1:dim
    if cxl(j)==0
        cxl(j)=0.1;
    end
    if cxl(j)==0.25
        cxl(j)=0.26;
    end
    if cxl(j)==0.5
        cxl(j)=0.51;
    end
    if cxl(j)==0.75
        cxl(j)=0.76;
    end
    if cxl(j)==1
        cxl(j)=0.9;
    end
end
for j=1:dim
    cxl(j)=4*cxl(j)*(1-cxl(j));        %logic混沌方程
end
% if Boundary_no==1
    Positions=cxl.*(ub-lb)+lb;
% end
% If each variable has a different lb and ub
% if Boundary_no>1
%     for i=1:dim
%         ub_i=ub(i);
%         lb_i=lb(i);
%         Positions(i,:)=cxl.*(ub_i-lb_i)+lb_i;
%     end
end