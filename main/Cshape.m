function newCenter = Cshape(I, oldCenter)
% ����mylinesegֱ�߷ָѡ��ͼ��������ֵ��߶��ϵ��е���Ϊ�µ�����
    if ~islogical(I)
        %������Ƕ�ֵͼ������ת��Ϊ��ֵͼ
        try
            I=rgb2gray(I);
        catch
        end
        level=graythresh(I);
        I=im2bw(I,level);
    end
    I0 = I;
%     theta0 = [1:6:180,179]; %������31��ֱ�ߣ�ѡ���ŵ�
    theta0 = 0:10:179; % �����ߵ����������Ч�� % 2015/4/11
    [m,n]=size(I); % m��n��
    Sr = 10000; % �����
    for ii = 1:length(theta0)
        I = I0;
        theta = theta0(ii);
        newCenter = zeros(n,2); %����ÿ��ֱ��ʱ����ֱ������Ľ������˾���
        % �����������ֱ�ߺ󣬴���ѡ���е㼴�ɡ�ÿ��xһ���㣬��๲n����
        for xx = 1:n
            % Ҫ��֤ÿһ��x���ж�Ӧ��y���ӣ��������ܵõ�����������
            besty = 1;
            deltay = 10000;
            tag = 0; %б����ӽ��ĸ����������������Ϊ1������Ϊ0��Ĭ�ϲ���������
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
                I(besty, xx) = 0; %����besty�ĸ���Ҫ��Ϊ0
            end
            if tag == 1
                % ����������ϵĵ㣬�����newCenter
                newCenter(xx,:) = [xx, besty];
            end
        end
        I(round(oldCenter(2)),round(oldCenter(1))) = 0;
        
        % ������ֱ�߷ָ���ͼ�������ж�����������
        se=strel('square',3); % ������С�ޣ�������С�ĽṹԪ�أ���СΪ3������2��Ч��û��3�á��������ʵ�����������
        I=imopen(I,se); % �ȸ�ʴ�����ͣ�����С���塢����ϸ���������塢ƽ��������߽�
% % %         I=bwareaopen(I,10,4); % ��ˮ��ָ����ַǳ���С�飬�˴�������ɺ�������(10fishring.avi frame=45)
        [L,num] = bwlabel(I, 4); %����ͨ����
        if num >= 2 
        %������Ϊnum==2ʱ�п��ܸ����з���û���ֳ�num=2��(10fishring.avi frame=45),�Ӷ����³���
            S0 = regionprops(L,'Area'); %����ͨ�������
            S0 = cat(1,S0.Area);
            S0 = sort(S0,'descend'); %�Ƚ������������
            Srtmp = max([S0(1)/S0(2), S0(2)/S0(1)]); % ���߼������
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