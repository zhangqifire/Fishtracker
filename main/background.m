function back=background(video, resize, numF)
% ʶ����Ƶ�еĹ̶������������ػҶȱ���ͼ����
% videoΪVideoReader����Ķ���, resizeΪ���ű���
    if nargin<2
        resize=0.5;
    end
    if nargin<3
        numF=30;        
    end
    numFrames=video.NumberOfFrames;
    p1=1; p2=numFrames;
    p=randi([p1,p2],numF,1);
    I=read(video,p(1));
    I=imresize(I,resize); %����
    I=rgb2gray(I);
    [row,col,dim]=size(I);
    II=zeros(row,col,dim,length(p));
    for j=2:length(p)
        I=read(video,p(j));
        I=imresize(I,resize); %����
        I=rgb2gray(I);
        II(:,:,:,j)=I;
    end
    back = median(double(II),4); %��λ��
%     back = mean(double(II),4); %��ֵ
%     back = mode(double(II),4); %����
    back = uint8(back);
%     imshow(back);
%     save('background_cm_0.4.mat','back');
%     imwrite(back,'background_cm_0.4.jpg','jpg');
end