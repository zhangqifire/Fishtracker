function [Imax, deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS)
    % ��Ծ������⣬�޸�ע�����С�����ʵ������������Ĳ��ֲ���
    % ���ص������ù��������ĵ��ֱ�߷ָѡ�ܽ�ͼ�ξ��ֵ�ֱ�ߣ��÷����Ƚ��������������ص������
    % ���ڳ����������ص��������׼ȷ���½�����С����ָ��
    if numtmp < 2.5
        deltaConf = 1;
    else
        deltaConf = 2;
    end
    Imax0 = Imax;
%     theta0 = [1:6:180,179]; %������31��ֱ�ߣ�ѡ���ŵ�
    theta0 = 0:10:179; % �����ߵ����������Ч��
    [m,n]=size(Imax); % m��n��
    Sr = 10000; % �����
    Ecr = 10000; % �����ʱ�
    Ecr1 = 10000; % �����ʱ�
    Ecr2 = 10000; % �����ʱ�
    Imaxtmp = Imax0; % ���ս��
    for ii = 1:length(theta0)
        Imax = Imax0;
        theta = theta0(ii);
        % 2015-3-17 ����ÿһ��yҲҪ�ж�Ӧ��x����
        % Ҫ��֤ÿһ��x���ж�Ӧ��y���ӣ�ͬʱÿһ����yҲҪ�ж�Ӧ��x���ӣ��������ܵõ�����������
        for yy = 1:m
            % ��֤ÿһ��y���ж�Ӧ��x����
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
            Imax(yy, bestx) = 0; %����bestx�ĸ���Ҫ��Ϊ0
        end
        for xx = 1:n
            % ��֤ÿһ��x���ж�Ӧ��y����
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
            Imax(besty, xx) = 0; %����besty�ĸ���Ҫ��Ϊ0
        end
        Imax(round(Cmax(2)),round(Cmax(1))) = 0;
%         figure;imshow(Imax)
        % ������ֱ�߷ָ���ͼ�������ж�������ֺ����������
        se=strel('square',3); % ������С�ޣ�������С�ĽṹԪ�أ���СΪ3������2��Ч��û��3�á��������ʵ�����������
        Imax=imopen(Imax,se); % �ȸ�ʴ�����ͣ�����С���塢����ϸ���������塢ƽ��������߽�
        Imax=bwareaopen(Imax,round(meanS/4),4);
        [L,num] = bwlabel(Imax, 4); %����ͨ����
        if num >= 2 
            % �˴�ԭλ num == 2,������ʱ����ͷ���ֳ��˶����,����δ������,�Ӷ������ѭ��
            % ��������Cshape��v3�з��ֵ�һ��,����ֱ��v4���޸�
            S0 = regionprops(L,'Area'); %����ͨ�������
            S0 = cat(1,S0.Area);
            [S0,indtmp] = sort(S0,'descend'); %�Ƚ������������
            Ectmp = regionprops(L,'Eccentricity'); %����ͨ�����������
            Ectmp = cat(1,Ectmp.Eccentricity);
            Ectmp = Ectmp(indtmp); % ���������������������
            Srtmp = max([S0(1)/S0(2), S0(2)/S0(1)]); % ���߼������
            Ecrtmp1 = max([Ectmp(1)/medianEc, medianEc/Ectmp(1)]); % ����λֵ�����ʱ�
            Ecrtmp2 = max([Ectmp(2)/medianEc, medianEc/Ectmp(2)]); % ����λֵ�����ʱ�
            Ecrtmp = max([Ectmp(1)/Ectmp(2), Ectmp(2)/Ectmp(1)]); % ���߼������ʱ�
            if Srtmp+(Ecrtmp+Ecrtmp1+Ecrtmp2)*5 < Sr+(Ecr+Ecr1+Ecr2)*5
                % �����ʵ�Ȩ��Ӧ�ô�Щ�����幫ʽ���Ը���ʵ���������
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