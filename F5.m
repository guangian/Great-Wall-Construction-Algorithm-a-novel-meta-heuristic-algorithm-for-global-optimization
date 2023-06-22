function [fit,result,x0]=F5(x)
x0=x;
global data
%% 解码获得厂区布局
x1=x(1:data.numUnit+data.maxS-1);
x(1:data.numUnit+data.maxS-1)=[];
X=x;
map=zeros(data.mapSize./data.accuracy);
[~,S]=sort(x1);
recording{1}=[];
len2=data.mapSize(2)+data.d(2);
len1=data.mapSize(1);
recording0.noUnit=S;
recording1=[];
jishu=0;

for i=1:length(S)
    no=S(i);
    if S(i)<=data.numUnit
        sizeUnit=data.Unit(no,1:2);
        type=1;
        jishu=jishu+1;
    else
        if jishu<data.minUnit
            continue;
        end
        sizeUnit=[data.r*2,data.r*2];
        type=2;
        jishu=0;
    end
    if ismember(no,[8,9]) && len1~=data.mapSize(1)
        len1=data.mapSize(1);
        len2=len2-max(recording{end}(:,5))-data.d(2);
        recording=[recording;{[]}];
        recording{end}=[recording{end};no,len1,len2,sizeUnit,type];
        recording1=[recording1;no,len1,len2,sizeUnit,type];
        len1=len1-(sizeUnit(1)+data.d(1));
        continue;
    end
    if ~isempty(recording1)
        position=find(recording1(:,1)==1 | recording1(:,1)==2);
        if ~isempty(position) && ismember(no,[1,2])
            len0=recording1(position,2);
            if len1<len0
                len1=len0;
                len2=len2-max(recording{end}(:,5))-data.d(2);
                recording=[recording;{[]}];
                recording{end}=[recording{end};no,len1,len2,sizeUnit,type];
                recording1=[recording1;no,len1,len2,sizeUnit,type];
                len1=len1-(sizeUnit(1)+data.d(1));
            else
                len1=len0;
                recording{end}=[recording{end};no,len1,len2,sizeUnit,type];
                recording1=[recording1;no,len1,len2,sizeUnit,type];
                len1=len1-(sizeUnit(1)+data.d(1));
            end
            continue;
        end
    end
    if len1-(sizeUnit(1)+data.d(1))>0
        recording{end}=[recording{end};no,len1,len2,sizeUnit,type];
        recording1=[recording1;no,len1,len2,sizeUnit,type];
        len1=len1-(sizeUnit(1)+data.d(1));
    else
        len1=data.mapSize(1);
        len2=len2-max(recording{end}(:,5))-data.d(2);
        recording=[recording;{[]}];
        recording1=[recording1;no,len1,len2,sizeUnit,type];
        recording{end}=[recording{end};no,len1,len2,sizeUnit,type];
        len1=len1-(sizeUnit(1)+data.d(1));
    end
end
%% 将各单元中心对齐
recording0.unit=[];
jishu1=1;
for i=1:length(recording)
    if rem(i,2)==0
        index=length(recording{i}(:,1)):-1:1;
        recording{i}=recording{i}(index,:);
    end
    maxY=max(recording{i}(:,5));
    meanY=mean(recording{i}(:,5));    
    for j=1:length(recording{i}(:,1))
%         if recording{i}(j,1)<=data.numUnit
%             no=recording{i}(j,1);
%             unitType=recording{i}(j,6);
%         else
%             no=recording{i}(j,1)-data.numUnit;
%             unitType=0;
%         end
        no=recording{i}(j,1);
        unitType=recording{i}(j,6);
        x=data.mapSize(1)-recording{i}(j,2);
        y=recording{i}(j,3)-data.d(2);
        unitX=recording{i}(j,4);
        unitY=recording{i}(j,5);
        
        recording{i}(j,7)=x;
        recording{i}(j,8)=y-maxY/2+unitY/2;
        y=recording{i}(j,8);
        if no>data.numUnit
            recording{i}(j,1)=jishu1;
            recording{i}(j,6)=2;
            jishu1=jishu1+1;
        end
    end 
    recording0.unit=[recording0.unit;[recording{i},ones(length(recording{i}(:,1)),1)*i]];
end
Zone=[];
index=1;
jishu=1;  %单元指针
jishu1=1; %区域编号
for i=1:length(recording0.unit(:,1))
    if recording0.unit(i,6)==2 
        temp={recording0.unit((index:i-1),1)};
        if ~isempty(temp{1})
            Zone=[Zone;jishu,jishu+length(temp{1})-1];
            index1=temp{1};
            Zone1(index1)=jishu1;
            Zone2{jishu1}=index1;
            jishu1=jishu1+1;
            jishu=jishu+length(temp{1});
        end
        index=i+1;
    end
end
temp={recording0.unit(index:end,1)};
if ~isempty(temp{1})
    Zone=[Zone;jishu,jishu+length(temp{1})-1];
    index1=temp{1};
    Zone1(index1)=jishu1;
    Zone2{jishu1}=index1;
end
position1=find(recording0.unit(:,6)==2);
if length(position1)>length(Zone(:,1))-1
    recording0.unit(position1(length(Zone(:,1)):end),:)=[];
end
position1=find(recording0.unit(:,6)==2);
position2=find(recording0.unit(:,6)==1);
recording0.unit=[recording0.unit(position2,:);recording0.unit(position1,:)];
recording0.noUnit=recording0.unit(:,1);
recording0.Zone=Zone;
recording0.Zone1=Zone1; %每个节点的区域
recording0.Zone2=Zone2; %每个区域的节点
% % 安排缓冲区
% for i=1:length(recording0.Zone(:,1))-1
%     no1=recording0.Zone(i,2);
%     no2=recording0.Zone(i+1,2);
%     layer1=recording0.unit(no1,9);
%     layer2=recording0.unit(no1,9);
%     x1=recording0.unit(no1,7);
%     x2=recording0.unit(no2,7);
%     y1=recording0.unit(no1,8);
%     y2=recording0.unit(no2,8);
%     unitX1=recording0.unit(no1,4);
%     unitX2=recording0.unit(no2,4);
%     unitY1=recording0.unit(no1,5);
%     unitY2=recording0.unit(no2,5);
%     if layer1==layer2
%         if rem(layer1,2)==1
%             x=(x1+unitX1+data.d(1)/2)-data.r;
%             y=y1-unitY1/2+data.r;
%         else
%             x=(x1-data.d(1)/2)-data.r;
%             y=y1-unitY1/2+data.r;
%         end
%     else
%         y=(y1+unitY1/2+y2-unitY2/2)/2+data.r;
%         x=x1;
%     end
%     temp=[i,0,0,data.r*2,data.r*2,2,x,y,0];
%     recording0.unit=[recording0.unit;temp];
% end
% 计算运载量和
data.numZone=max(recording0.Zone1);
DQ=zeros(1,data.numZone);
DQmax=zeros(1,data.numZone);
T=zeros(1,data.numZone);
for i=1:data.numUnit
    for j=1:data.numUnit
        Q=data.Distance(i,j);
        noZone1=recording0.Zone1(i);
        noZone2=recording0.Zone1(j);
        p1=find(recording0.unit(:,1)==i & recording0.unit(:,6)==1);
        p2=find(recording0.unit(:,1)==j & recording0.unit(:,6)==1);
        x1=recording0.unit(p1,7);
        y1=recording0.unit(p1,8);
        x2=recording0.unit(p2,7);
        y2=recording0.unit(p2,8);
        if noZone1==noZone2
            D=abs(x1-x2)+abs(y1-y2);
            temp=D*Q;
            DQ(noZone1)=DQ(noZone1)+temp;
            T(noZone1)=T(noZone1)+D/data.v+data.upLoadT+data.downLoadT;
        else
            if noZone1<noZone2
                index= noZone1:noZone2;
                biaoji=1;
            else
                index= noZone1:-1:noZone2;
                biaoji=2;
            end
            for k=1:length(index)
                if k==1
                    p1=find(recording0.unit(:,1)==i & recording0.unit(:,6)==1);
                    if biaoji==1
                        p2=find(recording0.unit(:,1)==index(k) & recording0.unit(:,6)==2);
                    else
                        p2=find(recording0.unit(:,1)==index(k+1) & recording0.unit(:,6)==2);
                    end
                elseif k==length(index)
                    p1=find(recording0.unit(:,1)==j & recording0.unit(:,6)==1);
                    if biaoji==1
                        p2=find(recording0.unit(:,1)==index(k-1) & recording0.unit(:,6)==2);
                    else
                        p2=find(recording0.unit(:,1)==index(k) & recording0.unit(:,6)==2);
                    end
                else
                    if biaoji==1
                        p1=find(recording0.unit(:,1)==index(k) & recording0.unit(:,6)==2);
                        p2=find(recording0.unit(:,1)==index(k-1) & recording0.unit(:,6)==2);
                    else
                        p1=find(recording0.unit(:,1)==index(k) & recording0.unit(:,6)==2);
                        p2=find(recording0.unit(:,1)==index(k+1) & recording0.unit(:,6)==2);
                    end
                end
                x1=recording0.unit(p1,7);
                y1=recording0.unit(p1,8);
                x2=recording0.unit(p2,7);
                y2=recording0.unit(p2,8);
                D=abs(x1-x2)+abs(y1-y2);
                temp=D*Q;
                DQ(index(k))=DQ(index(k))+temp;
                DQmax(index(k))=DQmax(index(k))+data.Dmax*Q;
                T(index(k))=T(index(k))+D/data.v+data.upLoadT+data.downLoadT; 
            end
        end
    end
end
if max(DQ)>data.maxDQ
    punishiment=max(DQ)-data.maxDQ;
else
    punishiment=0;
end
position=find(recording0.unit(:,8)-recording0.unit(:,5)<0);
punishiment=punishiment+1000*length(position);
if max(recording0.Zone1)==1
    punishiment=punishiment+1000000;
end
for i=1:max(recording0.Zone1)
    position=find(recording0.Zone1==i);
    if length( position)<data.minUnit
        punishiment=punishiment+1000000;
    end
end
numAGV=max(recording0.Zone1);
%%
minX=min(recording0.unit(:,7));
maxX=max(recording0.unit(:,7)+recording0.unit(:,4));
maxY=max(recording0.unit(:,8));
minY=min(recording0.unit(:,8)-recording0.unit(:,5));
fit1=sum(sum(DQ))/sum(sum(DQmax));
fit2=(maxX-minX)*(maxY-minY);
fit3=numAGV;
fit=sum(data.w.*[fit1,fit2/data.LWmax,fit3/data.maxAGV])+punishiment*10;
if nargout>1
    result.fit=fit;
    result.recording0=recording0;
    result.DQ=DQ;
end
end