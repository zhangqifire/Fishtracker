function [Imax, deltaConf] = mykmeans( Imax, nfish, tmpeps )
%MYKMEANS ��k��ֵ���෨�ָ�ͼ��
%   ����K-means���࣬�ҳ���֮��ķֽ��ߣ��ٰѷֽ��߱�ڣ����طָ���ͼ��
    if nargin < 3
        tmpeps=0.15; % ѡ��ֽ�����ֵ
    end
    L = bwlabel(Imax, 4); %����ͨ����
    Ls = sparse(L);
    [i,j]=find(Ls);
    [Idx, ~, ~, D] = kmeans([i,j],nfish);
    for ii=1:length(Idx)
        % ����ÿ���㣬�����2�����������ľ���ȽϽӽ�����˵�����ڽ��㴦��������
        minD = min(D(ii,:));
        % �����2�����ĵľ����С��2����С����С�����tmpeps�������ж�Ϊ�ֽ����ϵĵ�
        ind = find(D(ii,:)-minD < max(tmpeps * minD, 2)); 
        if length(ind)>1
            Imax(i(ii),j(ii)) = 0;
        end
    end    
    deltaConf = 3;
end

