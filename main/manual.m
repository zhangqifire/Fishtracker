%manual 人工分割连在一起的个体
%   通过在画面中用鼠标选择线段的两个端点，将线段上的所有点变成纯白（黑色个体）或纯黑（白色个体）
% 2015.6.11 减去背景后的图片总是鱼更亮（不管之前鱼是暗还是亮），所以分割线得是纯黑

if rgb == 1
    if strcmp(darkfish,'yes')
        tag = [255,255,255]; % for RGB,因为是神经网络识别，黑鱼，线应该白，白鱼线应该黑
    else
        tag = [0,0,0];
    end
else
    tag = 0; % for GRAY
end
prompt1 = {'go on?'};
name1 = 'Segmenting manual';
numlines1 = 1;
defaultanswer1={'yes'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
in = inputdlg(prompt1,name1,numlines1,defaultanswer1,options);
try
    go = in{1};
catch
    go = 'no';
end
while(strcmp(go,'yes')==1)
    imshow(I);
    title('choose two point for the segmentation and export cursor_info')
    cursor_info=1;
    try
        while(cursor_info==1)
            pause(1);
        end
    catch
    end
    pos=cat(1,cursor_info.Position);
    x1=pos(1,1);y1=pos(1,2);x2=pos(2,1);y2=pos(2,2);
    slope=(y2-y1)/(x2-x1);
    if slope==Inf
        I(x1,y1:y2,:)=tag;
    elseif slope==0
        I(x1:x2,y1,:)=tag;
    else
        if x1<x2
            tagx=1;
        else
            tagx=-1;
        end
        if y1<y2
            tagy=1;
        else
            tagy=-1;
        end
        % 下面保证至少每行有一个点
        for xx=x1:tagx:x2
            besty = y1;
            bestslope = Inf;
            for yy=y1:tagy:y2
                tmpslope = (yy-y1)/(xx-x1);
                if abs(tmpslope-slope) < abs(bestslope-slope)
                    bestslope = tmpslope;
                    besty = yy;
                end
                if abs(tmpslope-slope) < 0.02
                    I(yy,xx,:) = tag;
                end                
            end
            I(besty, xx,:) = tag; %至少besty的格子要变为0
        end
        % 下面保证每列至少有一个点，这样保证是一条连续的线
        for yy=y1:tagy:y2
            bestx = x1;
            bestslope = Inf;
            for xx=x1:tagx:x2
                tmpslope = (yy-y1)/(xx-x1);
                if abs(tmpslope-slope) < abs(bestslope-slope)
                    bestslope = tmpslope;
                    bestx = xx;
                end
                if abs(tmpslope-slope) < 0.02
                    I(yy,xx,:) = tag;
                end                
            end
            I(yy, bestx,:) = tag; %至少bestx的格子要变为0
        end
    end
    imshow(I);
    prompt1 = {'go on?'};
    name1 = 'Segmenting manual';
    numlines1 = 1;
    defaultanswer1={'yes'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    in = inputdlg(prompt1,name1,numlines1,defaultanswer1,options);
    try
        go = in{1};
    catch
        go = 'no';
    end
end
close(gcf)
