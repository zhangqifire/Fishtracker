function [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s2( I,C0indtmp,Conf0,num0,minS )
% ������Ķ�ֵͼ�е��ض���ͨ���Ž��зָ�

% 2015/9/7  ���ηָ��ͼ����Ȼ��С����С���һ�����ͨ����ȥ��
% 2015/5/1  �޸ķ�ˮ��ָ�Ϊk-means����ָ�
    
    % ����Եõ�������ͨ���Ž��зָ��ԭ��fishsegn_s1�����ָ�Ľ��
    [L,num] = bwlabel(I, 4); %����ͨ����
    S0 = regionprops(L,'Area'); %����ͨ�������
    S0 = cat(1,S0.Area);
    C0 = regionprops(L,'Centroid'); %����ͨ����λ��
    C0 = cat(1,C0.Centroid);
    len0 = regionprops(L,'MajorAxisLength'); % ��Բ����
    len0 = cat(1,len0.MajorAxisLength);
    wid0 = regionprops(L,'MinorAxisLength'); % ��Բ����
    wid0 = cat(1,wid0.MinorAxisLength);
    Ec0 = regionprops(L,'Eccentricity'); % ������
    Ec0 = cat(1,Ec0.Eccentricity);
    B0=regionprops(L,'BoundingBox'); % ����ͨ����߽�
    B0=cat(1,B0.BoundingBox);
%     O0=regionprops(L,'Orientation'); % ��Բ����
%     O0=cat(1,O0.Orientation);
    medianEc = median(Ec0); % ��������λ����ƽ��ֵ�ᱻ��ƫ
    meanS = sum(sum(I)) / num0; % ������ƽ��������ص���ʹ���ƫС������ʹ���ƫ��������ƫС�Ӷ�
    
    S0tmp = S0(C0indtmp);
    maxPos = C0indtmp(find(S0tmp==max(S0tmp),1,'first')); % ���������ͨ�����е����
    % ����Ը�����Ž��зָ�
    
% % %             len1 = len0(maxPos); % ��Բ����
% % %             wid1 = wid0(maxPos); % ��Բ����
    r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
    c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
    Imax = I(r1:r2,c1:c2); % �����ͨ���ž���ͼ������ָ�������־���
    numtmp = sum(sum(Imax))/meanS; % ���Ƶ��ص�������
% % %             lentmp = min(len0); widtmp = min(wid0);
    deltaConf = 0;
    if numtmp > 2.75
        % �ص�����϶࣬�÷�ˮ�뷨
        [I(r1:r2,c1:c2), deltaConf] = mykmeans(Imax,round(numtmp)); % K��ֵ����ָ��
%         [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % ��ˮ��ָ��
    elseif sum(sum(Imax)) > max(meanS, minS)
        % �������ٴ�����С����ŷָ�������ص��������С��minS�����Ժ�С������2֡���������Ŀ����Լ���Ϊ0
        Cmax = C0(maxPos,:) - B0(maxPos,1:2); % ��������ĵ�λ��,Ҫת��Ϊ��Imax��ͼ�е�����
        [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % ֱ�߷ָ��
    end
    Conf0 = Conf0 - deltaConf; % �ָ��ص��㽫����׼ȷ��
    
    % 2015.9.6  ���ηָ��ͼ����Ȼ��С����С���һ�����ͨ����ȥ��
    I=bwareaopen(I,max(1,floor(minS*0.5)),4); %��4�ھӲ�����ͨ���ţ�ȥ�����ص�����minS*0.5����
    
    % �ָ�����ͼ����������ͨ����
    [L,num] = bwlabel(I, 4); %����ͨ����
    S0 = regionprops(L,'Area'); %����ͨ�������
    S0 = cat(1,S0.Area);
    C0 = regionprops(L,'Centroid'); %����ͨ����λ��
    C0 = cat(1,C0.Centroid);
    B0=regionprops(L,'BoundingBox'); % ����ͨ����߽�
    B0=cat(1,B0.BoundingBox);           
    Ec0 = regionprops(L,'Eccentricity'); % ������
    Ec0 = cat(1,Ec0.Eccentricity);  
    len0 = regionprops(L,'MajorAxisLength'); % ��Բ����
    len0 = cat(1,len0.MajorAxisLength);
    wid0 = regionprops(L,'MinorAxisLength'); % ��Բ����
    wid0 = cat(1,wid0.MinorAxisLength);
%         O0=regionprops(L,'Orientation'); % ��Բ����
%         O0=cat(1,O0.Orientation);
    Iseg = I;
end