function [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s2( I,C0indtmp,Conf0,num0,minS )
% 对输入的二值图中的特定连通集团进行分割

% 2015/9/7  二次分割的图像依然把小于最小面积一半的联通集团去掉
% 2015/5/1  修改分水岭分割为k-means聚类分割
    
    % 下面对得到的鱼连通集团进行分割，还原出fishsegn_s1初步分割的结果
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
    
    S0tmp = S0(C0indtmp);
    maxPos = C0indtmp(find(S0tmp==max(S0tmp),1,'first')); % 最大团在连通集团中的序号
    % 下面对该最大团进行分割
    
% % %             len1 = len0(maxPos); % 椭圆长轴
% % %             wid1 = wid0(maxPos); % 椭圆短轴
    r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
    c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
    Imax = I(r1:r2,c1:c2); % 最大连通集团矩形图，下面分割这个部分就行
    numtmp = sum(sum(Imax))/meanS; % 估计的重叠鱼条数
% % %             lentmp = min(len0); widtmp = min(wid0);
    deltaConf = 0;
    if numtmp > 2.75
        % 重叠的鱼较多，用分水岭法
        [I(r1:r2,c1:c2), deltaConf] = mykmeans(Imax,round(numtmp)); % K均值聚类分割函数
%         [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % 分水岭分割函数
    elseif sum(sum(Imax)) > max(meanS, minS)
        % 必须至少大于最小面积才分割，两条鱼重叠后面积还小于minS可能性很小，连续2帧以上这样的可能性几乎为0
        Cmax = C0(maxPos,:) - B0(maxPos,1:2); % 最大团中心点位置,要转换为在Imax子图中的坐标
        [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % 直线分割函数
    end
    Conf0 = Conf0 - deltaConf; % 分割重叠鱼将降低准确度
    
    % 2015.9.6  二次分割的图像依然把小于最小面积一半的联通集团去掉
    I=bwareaopen(I,max(1,floor(minS*0.5)),4); %按4邻居查找连通集团，去除像素点少于minS*0.5个的
    
    % 分割后更新图像，重新找连通集团
    [L,num] = bwlabel(I, 4); %找连通区域
    S0 = regionprops(L,'Area'); %求连通区域面积
    S0 = cat(1,S0.Area);
    C0 = regionprops(L,'Centroid'); %求连通区域位置
    C0 = cat(1,C0.Centroid);
    B0=regionprops(L,'BoundingBox'); % 求连通区域边界
    B0=cat(1,B0.BoundingBox);           
    Ec0 = regionprops(L,'Eccentricity'); % 离心率
    Ec0 = cat(1,Ec0.Eccentricity);  
    len0 = regionprops(L,'MajorAxisLength'); % 椭圆长轴
    len0 = cat(1,len0.MajorAxisLength);
    wid0 = regionprops(L,'MinorAxisLength'); % 椭圆短轴
    wid0 = cat(1,wid0.MinorAxisLength);
%         O0=regionprops(L,'Orientation'); % 椭圆方向
%         O0=cat(1,O0.Orientation);
    Iseg = I;
end