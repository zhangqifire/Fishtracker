% function A = mycost(px1,py1,vx1,vy1,S1,ax1,ay1,len,px,py,S,v95,maxJump)
function [A,Av,Aa,As] = mycost(px1,py1,vx1,vy1,S1,ax1,ay1,len,px,py,S,v95,maxJump,countLost)
%MYCOST 根据输入信息计算成本矩阵
%   px1(t-1),py1(t-1),vx1(t-2),vy1(t-2),ax1(t-3),ay1(t-3),px(t),py(t)
%   要把[px(t),py(t)]与[px1(t-1),py1(t-1)]匹配
%   px1 M行1列，M为个体数目， px N行1列， N为识别出的个体数目， M与N可以不相等
%   通常 N >= M, 丢帧时 N < M.

%   2015/9/8  对于丢失的个体重新出现的，其可能在任何位置是对的，但是如果将成本设置为0则会导致在识别到的个体数量
%             少于实际个体数量时，丢失的个体占用未丢失个体的位置，而未丢失的个体变成新的丢失个体。
%             因此，A(ii,:)=9999999999(足够大的数，但不能是Inf)。
%             修改跳跃条件，连续丢帧数越大，可允许的搜索范围越大；
%             对于已经丢失的个体也要限定搜索范围，在可允许的搜索范围内的cost一样，之外的Inf
%   2015/4/11 修复权重，将面积变化的权重变小，只在速度和加速度差异小时才发生作用，因为面积容易因为误分割造成错误

    if nargin<12
        v95=40;
    end
    if nargin<13
        maxJump=1;
    end
    A = zeros(size(px1,1),size(px,1));
    
                %%% test %%%
    Av = zeros(size(px1,1),size(px,1));
    Aa = zeros(size(px1,1),size(px,1));
    As = zeros(size(px1,1),size(px,1));
                %%% test %%%
                
    wvel = 10; % 速度超出常规速度的权重，由此排除远的。
    wacce = 5; % 加速度的权重，由此排除连续交叉的，还应该考虑加速度鱼
    warea = 3; % 面积的权重，
    
%     vel = sqrt((vx1(:).^2 + vy1(:).^2)); % t-2时刻速度大小
    for ii = 1:size(px1,1)
        % 在鱼从未丢失到丢失及从丢失到再次找回时，与其相关的位置，速度，加速度都可能是NaN
        % 因此计算时需要单独处理，否则在调用C语言程序时会导致程序崩溃
        if isnan(px1(ii)) || isnan(vx1(ii)) || isnan(ax1(ii))
            % 如果个体在上一帧已经完全丢失了，则现在可以在任何位置
            A(ii,:) = 9999999999; % 所以其成本全部一样，且为9999999999
            % 2015.9.8 对于已经丢失的个体也要限定搜索范围，在可允许的搜索范围内的cost一样，之外的Inf
            for jj = 1:size(px,1)
                % 计算第ii个source与第jj个target匹配的成本，成本越小越好
                vx2 = px(jj) - px1(ii); 
                vy2 = py(jj) - py1(ii);
                v2 = sqrt(vx2^2+vy2^2);  %匹配后source的t-1时刻速度大小
                if v2 > v95 * (countLost(ii) + 2) + maxJump * len(ii) 
                    % 2015.9.8 修改跳跃条件，连续丢帧数越大，可允许的搜索范围越大
                    A(ii,jj) = A(ii,jj) + Inf;  % 如果鱼的速度太大，表明跳跃了
                end
            end
        else
            for jj = 1:size(px,1)
                % 计算第ii个source与第jj个target匹配的成本，成本越小越好
                vx2 = px(jj) - px1(ii); 
                vy2 = py(jj) - py1(ii);
                v2 = sqrt(vx2^2+vy2^2);  %匹配后source的t-1时刻速度大小
                ax2 = vx2 - vx1(ii); 
                ay2 = vy2 - vy1(ii);
                a2 = ax2^2 + ay2^2; % 匹配后source在t-2时刻的加速度
                a1 = ax1(ii)^2 + ay1(ii)^2; % 匹配后source在t-3时刻的加速度
                if v2 > v95
                    tmp = v2 - v95; % 速度超常越多，成本越大
                else
                    tmp = 0;
                end
                
                %%% test %%%
                Av(ii,jj) = tmp;
                Aa(ii,jj) = sqrt(a1 + a2);
                As(ii,jj) = sqrt(abs(S(jj)-S1(ii)));
                %%% test %%%
                
                A(ii,jj) = wvel * tmp + wacce * sqrt(a1 + a2) + warea * sqrt(abs(S(jj)-S1(ii)));
                if v2 > v95 * (countLost(ii) + 2) + maxJump * len(ii) 
                    % 2015.9.8 修改跳跃条件，连续丢帧数越大，可允许的搜索范围越大
                    A(ii,jj) = A(ii,jj) + Inf;  % 如果鱼的速度太大，表明跳跃了
                end
            end
        end
    end
end