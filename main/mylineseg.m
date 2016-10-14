function [Imax, deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS)
    % 针对具体问题，修改注释中有“根据实际情况调整”的部分参数
    % 将重叠的鱼用过现在中心点的直线分割，选能将图形均分的直线，该方法比较适用于两条鱼重叠的情况
    % 对于超过三条鱼重叠的情况，准确度下降，减小信心指数
    if numtmp < 2.5
        deltaConf = 1;
    else
        deltaConf = 2;
    end
    Imax0 = Imax;
%     theta0 = [1:6:180,179]; %考虑这31条直线，选最优的
    theta0 = 0:10:179; % 减少线的条数，提高效率
    [m,n]=size(Imax); % m行n列
    Sr = 10000; % 面积比
    Ecr = 10000; % 离心率比
    Ecr1 = 10000; % 离心率比
    Ecr2 = 10000; % 离心率比
    Imaxtmp = Imax0; % 最终结果
    for ii = 1:length(theta0)
        Imax = Imax0;
        theta = theta0(ii);
        % 2015-3-17 加入每一个y也要有对应的x格子
        % 要保证每一个x都有对应的y格子，同时每一个个y也要有对应的x格子，这样才能得到连续的线条
        for yy = 1:m
            % 保证每一个y都有对应的x格子
            bestx = 1;
            deltax = Inf;
            for xx = 1:n
                deltatmp = abs((yy - Cmax(2)) / (xx - Cmax(1)) - tan(theta/180*pi));
                if deltatmp < deltax
                    deltax = deltatmp;
                    bestx = xx;
                end
                if deltatmp < 0.02
                    Imax(yy,xx) = 0;
                end
            end
            Imax(yy, bestx) = 0; %至少bestx的格子要变为0
        end
        for xx = 1:n
            % 保证每一个x都有对应的y格子
            besty = 1;
            deltay = Inf;
            for yy = 1:m
                deltatmp = abs((yy - Cmax(2)) / (xx - Cmax(1)) - tan(theta/180*pi));
                if deltatmp < deltay
                    deltay = deltatmp;
                    besty = yy;
                end
                if deltatmp < 0.02
                    Imax(yy,xx) = 0;
                end
            end
            Imax(besty, xx) = 0; %至少besty的格子要变为0
        end
        Imax(round(Cmax(2)),round(Cmax(1))) = 0;
%         figure;imshow(Imax)
        % 以上用直线分割了图像，下面判断面积均分和离心率情况
        se=strel('square',3); % 下面是小修，所以用小的结构元素，大小为3。试了2，效果没有3好。具体根据实际情况调整。
        Imax=imopen(Imax,se); % 先腐蚀后膨胀，消除小物体、在纤细处分离物体、平滑大物体边界
        Imax=bwareaopen(Imax,round(meanS/4),4);
        [L,num] = bwlabel(Imax, 4); %找连通区域
        if num >= 2 
            % 此处原位 num == 2,但是有时候大块头被分成了多个块,导致未被处理,从而造成死循环
            % 这个问题跟Cshape在v3中发现的一样,但是直到v4才修复
            S0 = regionprops(L,'Area'); %求连通区域面积
            S0 = cat(1,S0.Area);
            [S0,indtmp] = sort(S0,'descend'); %比较面积最大的两个
            Ectmp = regionprops(L,'Eccentricity'); %新连通区域的离心率
            Ectmp = cat(1,Ectmp.Eccentricity);
            Ectmp = Ectmp(indtmp); % 按照面积排序结果重新排序
            Srtmp = max([S0(1)/S0(2), S0(2)/S0(1)]); % 两者间面积比
            Ecrtmp1 = max([Ectmp(1)/medianEc, medianEc/Ectmp(1)]); % 与中位值离心率比
            Ecrtmp2 = max([Ectmp(2)/medianEc, medianEc/Ectmp(2)]); % 与中位值离心率比
            Ecrtmp = max([Ectmp(1)/Ectmp(2), Ectmp(2)/Ectmp(1)]); % 两者间离心率比
            if Srtmp+(Ecrtmp+Ecrtmp1+Ecrtmp2)*5 < Sr+(Ecr+Ecr1+Ecr2)*5
                % 离心率的权重应该大些，具体公式可以根据实际情况调整
                Sr = Srtmp;
                Ecr1 = Ecrtmp1;
                Ecr2 = Ecrtmp2;
                Ecr = Ecrtmp;
                Imaxtmp = Imax;
            end
        end
    end
    Imax = Imaxtmp;    
end