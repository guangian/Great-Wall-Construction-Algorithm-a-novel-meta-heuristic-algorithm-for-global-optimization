function drawPc_F5(result1,option,data,str)
recording0=result1.recording0;
t = tiledlayout('flow');
t.TileSpacing = 'tight';
t.Padding = 'tight';
rectangle('Position',[0,0,data.mapSize(1),data.mapSize(2)]) ;
hold on
for i=1:length(recording0.unit(:,1))
    no=recording0.unit(i,1);
    x=recording0.unit(i,7);
    y=recording0.unit(i,8);
    unitX=recording0.unit(i,4);
    unitY=recording0.unit(i,5);
    unitType=recording0.unit(i,6);
    text(x+unitX/3,y-unitY/3,num2str(no))
    if unitType==1
        rectangle('Position',[x,y-unitY,unitX,unitY]) ;
    else
        rectangle('Position',[x,y-unitY,unitX,unitY],'Curvature',[1,1]) ;
    end
end
title([str,'目标:',num2str(result1.fit)])
end