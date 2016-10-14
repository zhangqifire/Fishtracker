function [Ec0,len0, wid0, S0, C0, Conf0, Iseg] = fishsegn_s1( I,num0,strelSize,minS,maxS,minLight )
% 根据灰度图像，背景图像和个体数目识别出每个个体位置和面积
% 本函数根据输入图像进行初步分割，只对面积明显过大的团块进行分割
%   I -- 输入灰度图像
%   back -- 背景
%   num0 -- 图中个体数目
%   len0 -- 椭圆拟合的长轴长度
%   wid0 -- 椭圆拟合的短轴长度
%   Ec0 -- 离心率
%   S0 -- 个体面积
%   C0 -- 个体形心
%   Conf0 -- 识别准确度的信心指数
%   minS -- 最小鱼的面积, maxS同理
%   darkfish -- 鱼比背景更亮为yes，否则为no
%   strelSize -- 形态学操作结构元素大小，一般可以取4

% 2015/9/8  将亮度小于15的点变为全黑，之后再用graythresh函数，效果很好
% 2015/9/7  修改灰度转二值图的level值为0.05
% 2015/5/1  修改分水岭分割为k-means聚类分割
% 2015/4/11 修复转为二值图错误及导致后面Imax分割错误的bug
% 2015/4/9 修改判定什么样的连通集团需要继续分割的方法
%          数量少时，面积大于当前画面均值1.5倍的都分割；数量足够时，只分割面积大于最大值maxS的
% 2015/4/8 重大修改，针对rgb情况，只留下num0,strelSize这两个参数
% 2015/3/16 不要手动设置阈值，参数中的grayThreshold无用了
    
    if nargin < 4
        minS = 20; % 视频分辨率不要太低，否则鱼与错误识别的部分傻傻分不清
    end
    if nargin < 5
        maxS = Inf; % 视频分辨率高没关系
    end
    if nargin < 6
        minLight = 30; % 把亮度低于minLight的点亮度改为0
    end
    
    % 转为二值图，做基本形态学处理
    % 转换的二值图可能非常不理想
    I(I<minLight) = 0; % 2015.9.8 将亮度小于minLight的点变为全黑，之后再用graythresh函数，效果很好
    level=graythresh(I);
    I=im2bw(I,level);
    
    if strelSize >= 1
        se=strel('square',strelSize); % strelSize应该根据鱼大小调整，一般可以为3或4
        %http://www.cnblogs.com/tornadomeet/archive/2012/03/20/2408086.html
        % 2015.6.11 蚂蚁用imclose，鱼用imopen
%         I=imopen(I,se); % 先腐蚀后膨胀，消除小物体、在纤细处分离物体、平滑大物体边界，此处主要是干掉了尾巴
        I=imclose(I,se); % 先膨胀后腐蚀，填充内部空洞、连接临近物体、平滑大物体边界，此处不适用
    end
    I=bwareaopen(I,minS,4); %按4邻居查找连通集团，去除像素点少于minS个的
    
    % 下面对得到的鱼连通集团进行正式分割
    [L,num] = bwlabel(I, 4); %找连通区域
    S0 = regionprops(L,'Area'); %求连通区域面积
    S0 = cat(1,S0.Area);
    C0 = regionprops(L,'Centroid'); %求连通区域位置
    C0 = cat(1,C0.Centroid);
    len0 = regionprops(L,'MajorAxisLength'); % 椭圆长轴
    len0 = cat(1,len0.MajorAxisLength);
    wid0 = regionprops(L,'MinorAxisLength'); % 椭圆短轴
    wid0 = cat(1,wid0.MinorAxisLength);
    Ec0 = regionprops(L,'Eccentricity'); % 离心率
    Ec0 = cat(1,Ec0.Eccentricity);
    B0=regionprops(L,'BoundingBox'); % 求连通区域边界
    B0=cat(1,B0.BoundingBox);
%     O0=regionprops(L,'Orientation'); % 椭圆方向
%     O0=cat(1,O0.Orientation);
    medianEc = median(Ec0); % 离心率中位数，平均值会被带偏
    meanS = sum(sum(I)) / num0; % 粗略鱼平均面积，重叠会使面积偏小，杂质使面积偏大，总体上偏小居多
    tagcut = 0; % 判定是否重新分割了
    if num == num0
        % 刚好每条鱼都是分离的,结果正确率最高
        Conf0 = 100; %信心指数为100
%         % 不分割了，分割导致块太多，减慢速度，降低精度，适得其反
%         if max(S0) > maxS
%             % 如果最大连通集团面积大于一定值，则把面积最大的几个连通集团分割
%             wait2cut = sum(S0>maxS);
%             while(wait2cut>0)
%                 tagcut = 1;
%                 wait2cut = wait2cut - 1;
%                 maxPos = find(S0==max(S0),1,'first'); % 最大团在连通集团中的序号
%     % % %             len1 = len0(maxPos); % 椭圆长轴
%     % % %             wid1 = wid0(maxPos); % 椭圆短轴
%                 r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
%                 c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
%                 Imax = I(r1:r2,c1:c2); % 最大连通集团矩形图，下面分割这个部分就行
%                 numtmp = sum(sum(Imax))/meanS; % 估计的重叠鱼条数
%     % % %             lentmp = min(len0); widtmp = min(wid0);
%                 if numtmp > 2.5
%                     % 重叠的鱼较多，用分水岭法
%                     [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % 分水岭分割函数
%                 else
%                     Cmax = C0(maxPos,:) - B0(maxPos,1:2); % 最大团中心点位置,要转换为在Imax子图中的坐标
%                     [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % 直线分割函数
%                 end
%                 Conf0 = Conf0 - deltaConf; % 分割重叠鱼将降低准确度
%                 % 分割后更新图像，重新找连通集团
%                 [L,num] = bwlabel(I, 4); %找连通区域
%                 S0 = regionprops(L,'Area'); %求连通区域面积
%                 S0 = cat(1,S0.Area);
%                 C0 = regionprops(L,'Centroid'); %求连通区域位置
%                 C0 = cat(1,C0.Centroid);
%                 B0=regionprops(L,'BoundingBox'); % 求连通区域边界
%                 B0=cat(1,B0.BoundingBox);   
%             end
%         end 
    elseif num > num0
        % 识别的鱼多于实际有的鱼，一般是有杂质没有去除，但后期连接各帧画面时一般不会出错
        Conf0 = 99; %信心指数
%         % 不分割了，分割导致块太多，减慢速度，降低精度，适得其反
%         if max(S0) > maxS
%             % 如果最大连通集团面积大于一定值，则把面积最大的几个连通集团分割
%             wait2cut = sum(S0>maxS);
%             while(wait2cut>0)
%                 tagcut = 1;
%                 wait2cut = wait2cut - 1;
%                 maxPos = find(S0==max(S0),1,'first'); % 最大团在连通集团中的序号
%     % % %             len1 = len0(maxPos); % 椭圆长轴
%     % % %             wid1 = wid0(maxPos); % 椭圆短轴
%                 r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
%                 c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
%                 Imax = I(r1:r2,c1:c2); % 最大连通集团矩形图，下面分割这个部分就行
%                 numtmp = sum(sum(Imax))/meanS; % 估计的重叠鱼条数
%     % % %             lentmp = min(len0); widtmp = min(wid0);
%                 if numtmp > 2.5
%                     % 重叠的鱼较多，用分水岭法
%                     [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % 分水岭分割函数
%                 else
%                     Cmax = C0(maxPos,:) - B0(maxPos,1:2); % 最大团中心点位置,要转换为在Imax子图中的坐标
%                     [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % 直线分割函数
%                 end
%                 Conf0 = Conf0 - deltaConf; % 分割重叠鱼将降低准确度
%                 % 分割后更新图像，重新找连通集团
%                 [L,num] = bwlabel(I, 4); %找连通区域
%                 S0 = regionprops(L,'Area'); %求连通区域面积
%                 S0 = cat(1,S0.Area);
%                 C0 = regionprops(L,'Centroid'); %求连通区域位置
%                 C0 = cat(1,C0.Centroid);
%                 B0=regionprops(L,'BoundingBox'); % 求连通区域边界
%                 B0=cat(1,B0.BoundingBox);   
%             end
%         end 
    else
        % 有鱼重叠或是未识别到
        Conf0 = 98;
        if num < num0 && max(S0) > maxS
            % 识别到的数目不足并且存在大于maxS集团，则把面积超规模的连通集团分割
%             wait2cut = num0 - num;
            wait2cut = sum(S0>maxS);
            while(wait2cut>0)
                tagcut = 1;
                wait2cut = wait2cut - 1;
                maxPos = find(S0==max(S0),1,'first'); % 最大团在连通集团中的序号
    % % %             len1 = len0(maxPos); % 椭圆长轴
    % % %             wid1 = wid0(maxPos); % 椭圆短轴
                
                r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
                c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
                % 2015/4/11,这里5Zebrafish_nocover_22min.avi_sf=1_ef=4725出错,索引超出矩阵维度。
                [r2Max,c2Max] = size(I);
                if r1<1
                    r1=1;
                end
                if r2>r2Max
                    r2=r2Max;
                end
                if c1<1
                    c1=1;
                end
                if c2>c2Max
                    c2=c2Max;
                end
                
                Imax = I(r1:r2,c1:c2); % 最大连通集团矩形图，下面分割这个部分就行
                numtmp = sum(sum(Imax))/meanS; % 估计的重叠鱼条数
    % % %             lentmp = min(len0); widtmp = min(wid0);
                if numtmp > 2.75
                    % 重叠的鱼较多，用分水岭法
                    [I(r1:r2,c1:c2), deltaConf] = mykmeans(Imax,round(numtmp)); % K均值聚类分割函数
%                     [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % 分水岭分割函数
                else
                    Cmax = C0(maxPos,:) - B0(maxPos,1:2); % 最大团中心点位置,要转换为在Imax子图中的坐标
                    [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % 直线分割函数
                end
                Conf0 = Conf0 - deltaConf; % 分割重叠鱼将降低准确度
                
                % 分割后更新图像，重新找连通集团
                [L,num] = bwlabel(I, 4); %找连通区域
                S0 = regionprops(L,'Area'); %求连通区域面积
                S0 = cat(1,S0.Area);
                C0 = regionprops(L,'Centroid'); %求连通区域位置
                C0 = cat(1,C0.Centroid);
                B0=regionprops(L,'BoundingBox'); % 求连通区域边界
                B0=cat(1,B0.BoundingBox);   
            end
        end 
    end
    
    % 如果重新分割了，就要重新统计
    if tagcut == 1 
        Ec0 = regionprops(L,'Eccentricity'); % 离心率
        Ec0 = cat(1,Ec0.Eccentricity);  
        len0 = regionprops(L,'MajorAxisLength'); % 椭圆长轴
        len0 = cat(1,len0.MajorAxisLength);
        wid0 = regionprops(L,'MinorAxisLength'); % 椭圆短轴
        wid0 = cat(1,wid0.MinorAxisLength);
%         O0=regionprops(L,'Orientation'); % 椭圆方向
%         O0=cat(1,O0.Orientation);
    end
    
    % 至此应该不存在未分割的大块头，但可能存在几乎重合的情况，所以鱼数量三种情况都可能
    for ii = 1:num
        % 对每条鱼判断是否中心点不在鱼身上，如果不在则用Cshape处理
        xx = round(C0(ii,1)); yy = round(C0(ii,2));
        if I(yy,xx)==0
            % 中心点不在鱼身上，用Cshape处理
            r1 = round(B0(ii,2)); r2 = round(B0(ii,2)+B0(ii,4));
            c1 = round(B0(ii,1)); c2 = round(B0(ii,1)+B0(ii,3));
            % 2015/4/11,这里5Zebrafish_nocover_22min.avi_sf=1_ef=4725出错,索引超出矩阵维度。
            [r2Max,c2Max] = size(I);
            if r1<1
                r1=1;
            end
            if r2>r2Max
                r2=r2Max;
            end
            if c1<1
                c1=1;
            end
            if c2>c2Max
                c2=c2Max;
            end
            Itmp = I(r1:r2,c1:c2); % C型鱼矩形图
            Ctmp = C0(ii,:) - B0(ii,1:2); % 小图中原始中心的坐标
%             Otmp = 90 - O0(ii); % 短轴与x轴的夹角，0-180度，一般是鱼身上的中心点在短轴方向上
            newC = Cshape(Itmp, Ctmp);
            C0(ii,:) = newC + B0(ii,1:2);
            Conf0 = Conf0 - 0.5; %中心点可能判断不准,但这种错误的影响相比多对多的影响小得多
        end
    end
    Iseg = I; %分割完成得到的二值图
end