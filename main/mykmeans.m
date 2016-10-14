function [Imax, deltaConf] = mykmeans( Imax, nfish, tmpeps )
%MYKMEANS 用k均值聚类法分割图像
%   先用K-means聚类，找出类之间的分界线，再把分界线变黑，返回分割后的图像
    if nargin < 3
        tmpeps=0.15; % 选择分界点的阈值
    end
    L = bwlabel(Imax, 4); %找连通区域
    Ls = sparse(L);
    [i,j]=find(Ls);
    [Idx, ~, ~, D] = kmeans([i,j],nfish);
    for ii=1:length(Idx)
        % 对于每个点，如果离2个及更多中心距离比较接近，则说明处于交汇处，将其变黑
        minD = min(D(ii,:));
        % 如果离2个中心的距离差小于2或者小于最小距离的tmpeps倍，就判定为分界线上的点
        ind = find(D(ii,:)-minD < max(tmpeps * minD, 2)); 
        if length(ind)>1
            Imax(i(ii),j(ii)) = 0;
        end
    end    
    deltaConf = 3;
end

