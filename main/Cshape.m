function newCenter = Cshape(I, oldCenter)
% 仿照mylineseg直线分割，选择将图形面积均分的线段上的中点作为新的中心
    if ~islogical(I)
        %如果不是二值图，则先转换为二值图
        try
            I=rgb2gray(I);
        catch
        end
        level=graythresh(I);
        I=im2bw(I,level);
    end
    I0 = I;
%     theta0 = [1:6:180,179]; %考虑这31条直线，选最优的
    theta0 = 0:10:179; % 减少线的条数，提高效率 % 2015/4/11
    [m,n]=size(I); % m行n列
    Sr = 10000; % 面积比
    for ii = 1:length(theta0)
        I = I0;
        theta = theta0(ii);
        newCenter = zeros(n,2); %尝试每条直线时，将直线与鱼的交点存入此矩阵，
        % 最后定下了最优直线后，从中选出中点即可。每个x一个点，最多共n个点
        for xx = 1:n
            % 要保证每一个x都有对应的y格子，这样才能得到连续的线条
            besty = 1;
            deltay = 10000;
            tag = 0; %斜率最接近的格子如果是鱼身上则为1，否则为0；默认不在鱼身上
            for yy = 1:m
                deltatmp = abs((yy - oldCenter(2)) / (xx - oldCenter(1)) - tan(theta/180*pi));
                if deltatmp < deltay
                    deltay = deltatmp;
                    besty = yy;
                end
                if deltatmp < 0.02
                    if I(yy,xx) == 1
                        tag = 1;
                        I(yy,xx) = 0;
                    end
                end
            end
            if I(besty, xx) == 1
                tag = 1;
                I(besty, xx) = 0; %至少besty的格子要变为0
            end
            if tag == 1
                % 如果是鱼身上的点，则存入newCenter
                newCenter(xx,:) = [xx, besty];
            end
        end
        I(round(oldCenter(2)),round(oldCenter(1))) = 0;
        
        % 以上用直线分割了图像，下面判断面积均分情况
        se=strel('square',3); % 下面是小修，所以用小的结构元素，大小为3。试了2，效果没有3好。具体根据实际情况调整。
        I=imopen(I,se); % 先腐蚀后膨胀，消除小物体、在纤细处分离物体、平滑大物体边界
% % %         I=bwareaopen(I,10,4); % 分水岭分割后出现非常多小块，此处容易造成后续出错(10fishring.avi frame=45)
        [L,num] = bwlabel(I, 4); %找连通区域
        if num >= 2 
        %当条件为num==2时有可能各种切法都没法分成num=2，(10fishring.avi frame=45),从而导致出错
            S0 = regionprops(L,'Area'); %求连通区域面积
            S0 = cat(1,S0.Area);
            S0 = sort(S0,'descend'); %比较面积最大的两个
            Srtmp = max([S0(1)/S0(2), S0(2)/S0(1)]); % 两者间面积比
            if Srtmp < Sr
                Sr = Srtmp;
                Centertmp = newCenter(newCenter(:,1)~=0,:);
            end
        end
    end
    try
        newCenter = Centertmp(ceil(length(Centertmp)/2),:);
    catch
        newCenter = oldCenter;
    end
end