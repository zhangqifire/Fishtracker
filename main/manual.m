%manual �˹��ָ�����һ��ĸ���
%   ͨ���ڻ����������ѡ���߶ε������˵㣬���߶��ϵ����е��ɴ��ף���ɫ���壩�򴿺ڣ���ɫ���壩
% 2015.6.11 ��ȥ�������ͼƬ���������������֮ǰ���ǰ��������������Էָ��ߵ��Ǵ���

if rgb == 1
    if strcmp(darkfish,'yes')
        tag = [255,255,255]; % for RGB,��Ϊ��������ʶ�𣬺��㣬��Ӧ�ðף�������Ӧ�ú�
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
        % ���汣֤����ÿ����һ����
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
            I(besty, xx,:) = tag; %����besty�ĸ���Ҫ��Ϊ0
        end
        % ���汣֤ÿ��������һ���㣬������֤��һ����������
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
            I(yy, bestx,:) = tag; %����bestx�ĸ���Ҫ��Ϊ0
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
