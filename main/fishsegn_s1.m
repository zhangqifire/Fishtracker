function [Ec0,len0, wid0, S0, C0, Conf0, Iseg] = fishsegn_s1( I,num0,strelSize,minS,maxS,minLight )
% ���ݻҶ�ͼ�񣬱���ͼ��͸�����Ŀʶ���ÿ������λ�ú����
% ��������������ͼ����г����ָֻ��������Թ�����ſ���зָ�
%   I -- ����Ҷ�ͼ��
%   back -- ����
%   num0 -- ͼ�и�����Ŀ
%   len0 -- ��Բ��ϵĳ��᳤��
%   wid0 -- ��Բ��ϵĶ��᳤��
%   Ec0 -- ������
%   S0 -- �������
%   C0 -- ��������
%   Conf0 -- ʶ��׼ȷ�ȵ�����ָ��
%   minS -- ��С������, maxSͬ��
%   darkfish -- ��ȱ�������Ϊyes������Ϊno
%   strelSize -- ��̬ѧ�����ṹԪ�ش�С��һ�����ȡ4

% 2015/9/8  ������С��15�ĵ��Ϊȫ�ڣ�֮������graythresh������Ч���ܺ�
% 2015/9/7  �޸ĻҶ�ת��ֵͼ��levelֵΪ0.05
% 2015/5/1  �޸ķ�ˮ��ָ�Ϊk-means����ָ�
% 2015/4/11 �޸�תΪ��ֵͼ���󼰵��º���Imax�ָ�����bug
% 2015/4/9 �޸��ж�ʲô������ͨ������Ҫ�����ָ�ķ���
%          ������ʱ��������ڵ�ǰ�����ֵ1.5���Ķ��ָ�����㹻ʱ��ֻ�ָ�����������ֵmaxS��
% 2015/4/8 �ش��޸ģ����rgb�����ֻ����num0,strelSize����������
% 2015/3/16 ��Ҫ�ֶ�������ֵ�������е�grayThreshold������
    
    if nargin < 4
        minS = 20; % ��Ƶ�ֱ��ʲ�Ҫ̫�ͣ������������ʶ��Ĳ���ɵɵ�ֲ���
    end
    if nargin < 5
        maxS = Inf; % ��Ƶ�ֱ��ʸ�û��ϵ
    end
    if nargin < 6
        minLight = 30; % �����ȵ���minLight�ĵ����ȸ�Ϊ0
    end
    
    % תΪ��ֵͼ����������̬ѧ����
    % ת���Ķ�ֵͼ���ܷǳ�������
    I(I<minLight) = 0; % 2015.9.8 ������С��minLight�ĵ��Ϊȫ�ڣ�֮������graythresh������Ч���ܺ�
    level=graythresh(I);
    I=im2bw(I,level);
    
    if strelSize >= 1
        se=strel('square',strelSize); % strelSizeӦ�ø������С������һ�����Ϊ3��4
        %http://www.cnblogs.com/tornadomeet/archive/2012/03/20/2408086.html
        % 2015.6.11 ������imclose������imopen
%         I=imopen(I,se); % �ȸ�ʴ�����ͣ�����С���塢����ϸ���������塢ƽ��������߽磬�˴���Ҫ�Ǹɵ���β��
        I=imclose(I,se); % �����ͺ�ʴ������ڲ��ն��������ٽ����塢ƽ��������߽磬�˴�������
    end
    I=bwareaopen(I,minS,4); %��4�ھӲ�����ͨ���ţ�ȥ�����ص�����minS����
    
    % ����Եõ�������ͨ���Ž�����ʽ�ָ�
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
    tagcut = 0; % �ж��Ƿ����·ָ���
    if num == num0
        % �պ�ÿ���㶼�Ƿ����,�����ȷ�����
        Conf0 = 100; %����ָ��Ϊ100
%         % ���ָ��ˣ��ָ�¿�̫�࣬�����ٶȣ����;��ȣ��ʵ��䷴
%         if max(S0) > maxS
%             % ��������ͨ�����������һ��ֵ�����������ļ�����ͨ���ŷָ�
%             wait2cut = sum(S0>maxS);
%             while(wait2cut>0)
%                 tagcut = 1;
%                 wait2cut = wait2cut - 1;
%                 maxPos = find(S0==max(S0),1,'first'); % ���������ͨ�����е����
%     % % %             len1 = len0(maxPos); % ��Բ����
%     % % %             wid1 = wid0(maxPos); % ��Բ����
%                 r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
%                 c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
%                 Imax = I(r1:r2,c1:c2); % �����ͨ���ž���ͼ������ָ�������־���
%                 numtmp = sum(sum(Imax))/meanS; % ���Ƶ��ص�������
%     % % %             lentmp = min(len0); widtmp = min(wid0);
%                 if numtmp > 2.5
%                     % �ص�����϶࣬�÷�ˮ�뷨
%                     [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % ��ˮ��ָ��
%                 else
%                     Cmax = C0(maxPos,:) - B0(maxPos,1:2); % ��������ĵ�λ��,Ҫת��Ϊ��Imax��ͼ�е�����
%                     [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % ֱ�߷ָ��
%                 end
%                 Conf0 = Conf0 - deltaConf; % �ָ��ص��㽫����׼ȷ��
%                 % �ָ�����ͼ����������ͨ����
%                 [L,num] = bwlabel(I, 4); %����ͨ����
%                 S0 = regionprops(L,'Area'); %����ͨ�������
%                 S0 = cat(1,S0.Area);
%                 C0 = regionprops(L,'Centroid'); %����ͨ����λ��
%                 C0 = cat(1,C0.Centroid);
%                 B0=regionprops(L,'BoundingBox'); % ����ͨ����߽�
%                 B0=cat(1,B0.BoundingBox);   
%             end
%         end 
    elseif num > num0
        % ʶ��������ʵ���е��㣬һ����������û��ȥ�������������Ӹ�֡����ʱһ�㲻�����
        Conf0 = 99; %����ָ��
%         % ���ָ��ˣ��ָ�¿�̫�࣬�����ٶȣ����;��ȣ��ʵ��䷴
%         if max(S0) > maxS
%             % ��������ͨ�����������һ��ֵ�����������ļ�����ͨ���ŷָ�
%             wait2cut = sum(S0>maxS);
%             while(wait2cut>0)
%                 tagcut = 1;
%                 wait2cut = wait2cut - 1;
%                 maxPos = find(S0==max(S0),1,'first'); % ���������ͨ�����е����
%     % % %             len1 = len0(maxPos); % ��Բ����
%     % % %             wid1 = wid0(maxPos); % ��Բ����
%                 r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
%                 c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
%                 Imax = I(r1:r2,c1:c2); % �����ͨ���ž���ͼ������ָ�������־���
%                 numtmp = sum(sum(Imax))/meanS; % ���Ƶ��ص�������
%     % % %             lentmp = min(len0); widtmp = min(wid0);
%                 if numtmp > 2.5
%                     % �ص�����϶࣬�÷�ˮ�뷨
%                     [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % ��ˮ��ָ��
%                 else
%                     Cmax = C0(maxPos,:) - B0(maxPos,1:2); % ��������ĵ�λ��,Ҫת��Ϊ��Imax��ͼ�е�����
%                     [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % ֱ�߷ָ��
%                 end
%                 Conf0 = Conf0 - deltaConf; % �ָ��ص��㽫����׼ȷ��
%                 % �ָ�����ͼ����������ͨ����
%                 [L,num] = bwlabel(I, 4); %����ͨ����
%                 S0 = regionprops(L,'Area'); %����ͨ�������
%                 S0 = cat(1,S0.Area);
%                 C0 = regionprops(L,'Centroid'); %����ͨ����λ��
%                 C0 = cat(1,C0.Centroid);
%                 B0=regionprops(L,'BoundingBox'); % ����ͨ����߽�
%                 B0=cat(1,B0.BoundingBox);   
%             end
%         end 
    else
        % �����ص�����δʶ��
        Conf0 = 98;
        if num < num0 && max(S0) > maxS
            % ʶ�𵽵���Ŀ���㲢�Ҵ��ڴ���maxS���ţ�����������ģ����ͨ���ŷָ�
%             wait2cut = num0 - num;
            wait2cut = sum(S0>maxS);
            while(wait2cut>0)
                tagcut = 1;
                wait2cut = wait2cut - 1;
                maxPos = find(S0==max(S0),1,'first'); % ���������ͨ�����е����
    % % %             len1 = len0(maxPos); % ��Բ����
    % % %             wid1 = wid0(maxPos); % ��Բ����
                
                r1 = round(B0(maxPos,2)); r2 = round(B0(maxPos,2)+B0(maxPos,4));
                c1 = round(B0(maxPos,1)); c2 = round(B0(maxPos,1)+B0(maxPos,3));
                % 2015/4/11,����5Zebrafish_nocover_22min.avi_sf=1_ef=4725����,������������ά�ȡ�
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
                
                Imax = I(r1:r2,c1:c2); % �����ͨ���ž���ͼ������ָ�������־���
                numtmp = sum(sum(Imax))/meanS; % ���Ƶ��ص�������
    % % %             lentmp = min(len0); widtmp = min(wid0);
                if numtmp > 2.75
                    % �ص�����϶࣬�÷�ˮ�뷨
                    [I(r1:r2,c1:c2), deltaConf] = mykmeans(Imax,round(numtmp)); % K��ֵ����ָ��
%                     [I(r1:r2,c1:c2), deltaConf] = mywatershed(Imax); % ��ˮ��ָ��
                else
                    Cmax = C0(maxPos,:) - B0(maxPos,1:2); % ��������ĵ�λ��,Ҫת��Ϊ��Imax��ͼ�е�����
                    [I(r1:r2,c1:c2), deltaConf] = mylineseg(Imax,numtmp,Cmax,medianEc,meanS); % ֱ�߷ָ��
                end
                Conf0 = Conf0 - deltaConf; % �ָ��ص��㽫����׼ȷ��
                
                % �ָ�����ͼ����������ͨ����
                [L,num] = bwlabel(I, 4); %����ͨ����
                S0 = regionprops(L,'Area'); %����ͨ�������
                S0 = cat(1,S0.Area);
                C0 = regionprops(L,'Centroid'); %����ͨ����λ��
                C0 = cat(1,C0.Centroid);
                B0=regionprops(L,'BoundingBox'); % ����ͨ����߽�
                B0=cat(1,B0.BoundingBox);   
            end
        end 
    end
    
    % ������·ָ��ˣ���Ҫ����ͳ��
    if tagcut == 1 
        Ec0 = regionprops(L,'Eccentricity'); % ������
        Ec0 = cat(1,Ec0.Eccentricity);  
        len0 = regionprops(L,'MajorAxisLength'); % ��Բ����
        len0 = cat(1,len0.MajorAxisLength);
        wid0 = regionprops(L,'MinorAxisLength'); % ��Բ����
        wid0 = cat(1,wid0.MinorAxisLength);
%         O0=regionprops(L,'Orientation'); % ��Բ����
%         O0=cat(1,O0.Orientation);
    end
    
    % ����Ӧ�ò�����δ�ָ�Ĵ��ͷ�������ܴ��ڼ����غϵ�����������������������������
    for ii = 1:num
        % ��ÿ�����ж��Ƿ����ĵ㲻�������ϣ������������Cshape����
        xx = round(C0(ii,1)); yy = round(C0(ii,2));
        if I(yy,xx)==0
            % ���ĵ㲻�������ϣ���Cshape����
            r1 = round(B0(ii,2)); r2 = round(B0(ii,2)+B0(ii,4));
            c1 = round(B0(ii,1)); c2 = round(B0(ii,1)+B0(ii,3));
            % 2015/4/11,����5Zebrafish_nocover_22min.avi_sf=1_ef=4725����,������������ά�ȡ�
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
            Itmp = I(r1:r2,c1:c2); % C�������ͼ
            Ctmp = C0(ii,:) - B0(ii,1:2); % Сͼ��ԭʼ���ĵ�����
%             Otmp = 90 - O0(ii); % ������x��ļнǣ�0-180�ȣ�һ���������ϵ����ĵ��ڶ��᷽����
            newC = Cshape(Itmp, Ctmp);
            C0(ii,:) = newC + B0(ii,1:2);
            Conf0 = Conf0 - 0.5; %���ĵ�����жϲ�׼,�����ִ����Ӱ����ȶ�Զ��Ӱ��С�ö�
        end
    end
    Iseg = I; %�ָ���ɵõ��Ķ�ֵͼ
end