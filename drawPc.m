function drawPc(result1,option,data,str,BestData)
figure
t = tiledlayout(1,2);
t.TileSpacing = 'tight';
t.Padding = 'tight';
% nexttile
% plot(BestData.Congervence,'LineWidth',2)
legend('PGWCA')
title('fitness curve')
nexttile
plot3(data.S0(:,1)*data.unit(1),data.S0(:,2)*data.unit(2),data.S0(:,3)*data.unit(3),'o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
hold on
plot3(data.E0(:,1)*data.unit(1),data.E0(:,2)*data.unit(2),data.E0(:,3)*data.unit(3),'h','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
plot3(result1.path(:,1).*data.unit(1),result1.path(:,2).*data.unit(2),result1.path(:,3).*data.unit(3),'-','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
for i=1:data.numObstacles
    x=1+data.Obstacle(i,1);
    y=1+data.Obstacle(i,2);
    z=1+data.Obstacle(i,3);
    long=data.Obstacle(i,4);
    wide=data.Obstacle(i,5);
    pretty=data.Obstacle(i,6);
    
    x0=ceil(x/data.unit(1))*data.unit(1);
    y0=ceil(y/data.unit(2))*data.unit(2);
    z0=ceil(z/data.unit(3))*data.unit(3);
    long0=ceil(long/data.unit(1))*data.unit(1);
    wide0=ceil(wide/data.unit(2))*data.unit(2);
    pretty0=ceil(pretty/data.unit(3))*data.unit(3);
    [V,F] = DrawCuboid(long0, wide0, pretty0, x0,y0,z0);
end
legend('起点','终点')
grid on
%axis equal
xlabel('x（km）')
ylabel('y（km）')
zlabel('z（km）')
title([str,'结果,最优目标:',num2str(result1.fit)])
% figure
% plot3(data.S0(:,1)*data.unit(1),data.S0(:,2)*data.unit(2),data.S0(:,3)*data.unit(3),'o','LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','r',...
%     'MarkerSize',10)
% hold on
% plot3(data.E0(:,1)*data.unit(1),data.E0(:,2)*data.unit(2),data.E0(:,3)*data.unit(3),'h','LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','r',...
%     'MarkerSize',10)
% plot3(result1.path(:,1).*data.unit(1),result1.path(:,2).*data.unit(2),result1.path(:,3).*data.unit(3),'-','LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','r',...
%     'MarkerSize',10)
% for i=1:data.numObstacles
%     x=1+data.Obstacle(i,1);
%     y=1+data.Obstacle(i,2);
%     z=1+data.Obstacle(i,3);
%     long=data.Obstacle(i,4);
%     wide=data.Obstacle(i,5);
%     pretty=data.Obstacle(i,6);
%     
%     x0=ceil(x/data.unit(1))*data.unit(1);
%     y0=ceil(y/data.unit(2))*data.unit(2);
%     z0=ceil(z/data.unit(3))*data.unit(3);
%     long0=ceil(long/data.unit(1))*data.unit(1);
%     wide0=ceil(wide/data.unit(2))*data.unit(2);
%     pretty0=ceil(pretty/data.unit(3))*data.unit(3);
%     [V,F] = DrawCuboid(long0, wide0, pretty0, x0,y0,z0);
% end
% legend('起点','终点')
% grid on
% xlabel('x（km）')
% ylabel('y（km）')
% zlabel('z（km）')
% title([str,'结果,最优目标:',num2str(result1.fit)])
end