function [fit,result,x0]=F4(x)
global data
    x0=x;
    %% 解码 获取每一个车辆的路径
    [~,S]=sort(x);
    path=cell(data.maxV,1);
    index=1;
    jishu=1;
    for i=1:length(x)
        if S(i)>data.numNode-1
            if ~isempty(S(index:i-1))
                path{jishu}=[1,S(index:i-1)+1,1];
            end
            index=i+1;
            jishu=jishu+1;
        end
    end
    if ~isempty(S(index:end))
        path{jishu}=[1,S(index:end)+1,1];
    end
    %%
    fit=0;
    for i=1:length(path)
        Load=0;
        recording.Path{i}=[];
        if ~isempty(path{i})
            t=0;
            Load=sum(data.node(path{i},4))+sum(data.node(path{i},5));
            for j=1:length(path{i})-1
                no1=path{i}(j);
                no2=path{i}(j+1);
                d=data.D(no1,no2);
                t0=d/data.v*60;
                w_M=data.node(no2,4);
                w_MB=data.node(no2,5);
                w_RMB=data.node(no2,6);
                num_M=w_M/data.W_milk;
                num_MB=w_MB/data.W_milk_B;
                num_RMB=w_RMB/data.W_milk_B;
                position=find(data.T1(:,1)<=num_M & data.T1(:,2)>=num_M);
                if ~isempty(position)
                    t1=data.T1(position,3);
                else
                    t1=0;
                end
                position=find(data.T2(:,1)<=num_MB & data.T2(:,2)>=num_MB);
                if ~isempty(position)
                    t2=data.T2(position,3);
                else
                    t2=0;
                end
                
                position=find(data.T2(:,1)<=num_RMB & data.T2(:,2)>=num_RMB);
                if ~isempty(position)
                    t3=data.T2(position,3);
                else
                    t3=0;
                end
                tt=t+t0+t1+t2+t3;
                C1=data.Cz*d;
                C2=data.Cw*d;
                C3=data.Cr*(t0+t1+t2+t3);
                C4=data.Cy*(data.Qk1+data.Qb1*Load);
                Cy=data.sigma*data.Qy*t0;
                Cq2=data.sigma*data.Qz*(t1+t2+t3);
                C5=Cy+Cq2;
                C6=data.P*data.b;
                C=C1+C2+C3+C4+C5+C6;
                recording.Path{i}=[recording.Path{i};no1,no2,d,t0,t1,t2,t3,t,tt,Load,C];
                fit=fit+C;
                t=tt;
                Load=Load-w_M-w_MB+w_RMB;
            end
            recording.T(i)=t;
            recording.Load(i)=max(recording.Path{i}(:,10));
            recording.numCustomer(i)=length(path{i})-2;
        else
            recording.T(i)=0;
            recording.Load(i)=0;
            recording.numCustomer(i)=0;
        end
    end
    %% 检查约束
    temp=recording.T-data.maxT;
    temp(temp<0)=0;
    punishiment1=sum(temp);
    temp=recording.Load-data.Q1;
    temp(temp<0)=0;
    punishiment2=sum(temp);
    temp=recording.numCustomer-data.maxCustom;
    temp(temp<0)=0;
    punishiment3=sum(temp);
    fit=fit+10000*(punishiment1+punishiment2+punishiment3);
    if nargout>1
        result.fit=fit;
        result.recording=recording;
    end
end