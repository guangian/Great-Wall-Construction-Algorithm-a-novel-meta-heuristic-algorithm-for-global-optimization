function drawPc_F4(result1,option,data,str)
    figure
    hold on
    legendStr=[{'车场'},{'顾客'}];
    plot(data.node(data.noCenter,2),data.node(data.noCenter,3),'h','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10);
    plot(data.node(data.noNode,2),data.node(data.noNode,3),'o','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10);
    for i=1:length(result1.recording.Path)
        path=[result1.recording.Path{i}(:,1);1];
        plot(data.node(path,2),data.node(path,3),'-','LineWidth',2);
        legendStr=[legendStr,{['第',num2str(i),'辆车路线']}];
    end
    legend(legendStr);
    title([str,'，求解路线，总目标：',num2str(result1.fit)]);
    for i=1:length(result1.recording.Path)
        figure
        hold on
        legendStr=[{'车场'},{'顾客'}];
        plot(data.node(data.noCenter,2),data.node(data.noCenter,3),'h','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','r',...
            'MarkerSize',10);
        plot(data.node(data.noNode,2),data.node(data.noNode,3),'o','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
        path=[result1.recording.Path{i}(:,1);1];
        plot(data.node(path,2),data.node(path,3),'-','LineWidth',2);
        legendStr=[legendStr,{['第',num2str(i),'辆车路线']}];
        legend(legendStr);
        title([str,'，第',num2str(i),'辆车路线，总目标：',num2str(result1.fit)]);
    end

end