%% Begin ANN Test
clc, close all
neighbors=8;
resize = 1;
[ffile,fpath] = uigetfile({'*.mp4';'*.avi'},'Select the Video file');
fname = [fpath,ffile];
% ��ȡͼ���ɰ�
region=imread([fname,'_�ɰ�.png']);
region=imresize(region,resize); %����

try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % ��ȡ��Ƶǰ���ж���Ƶ�Ƿ��Ѿ�����ȡ�������ڵ���ʱ���ٵȴ���ȡ��Ƶ��ʱ��
    end
catch
    V=VideoReader(fname); 
end

for kk=1:10
    frame = randi(V.NumberOfFrames);
    I0=read(V,frame); % ԭʼ��Ƶ֡ͼ��
    I0=imresize(I0,resize); %����   
    I=I0;
    % �ⲿ���Ǹ��ݺ�ɫ���ھ���ɫ�ݶ���Ϣ��ȡ���λ��
    % 2015.4.29 �����ھ���Ϣ����������ͼ�ξ�������bpfish.c�����������㷨�ж�
    Itmp = zeros(size(I,1)+4,size(I,2)+4,size(I,3));
    Itmp(3:end-2,3:end-2,:)=double(I);
    switch neighbors
        case 0
            rgbtmp = I;
            rgbtmp = permute(rgbtmp,[3 2 1]);
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2)*size(rgbtmp,3));
            % ��ʱrgbtmpÿһ�ж�ӦI��һ������������㷨���룬I(m,n)��Ӧ��n+(m-1)*size(I,2)
            I = myAnn0(rgbtmp);
        case 8
            r = double(I(:,:,1)); 
            rgbtmp = zeros(size(I,1),size(I,2),size(I,3),9);
            for jj=0:8
                rgbtmp(:,:,:,jj+1) = Itmp(2+floor(jj/3):end-3+floor(jj/3),2+mod(jj,3):end-3+mod(jj,3),:);
            end 
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2),size(rgbtmp,3)*size(rgbtmp,4));
            rgbtmp = permute(rgbtmp,[3 2 1]);
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2)*size(rgbtmp,3));
            % ��ʱrgbtmpÿһ�ж�ӦI��һ������������㷨���룬I(m,n)��Ӧ��n+(m-1)*size(I,2)
            I = myAnn8(rgbtmp);
        case 24
            r = double(I(:,:,1)); 
            rgbtmp = zeros(size(I,1),size(I,2),size(I,3),25);
            for jj=0:24
                rgbtmp(:,:,:,jj+1) = Itmp(1+floor(jj/5):end-4+floor(jj/5),1+mod(jj,5):end-4+mod(jj,5),:);
            end   
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2),size(rgbtmp,3)*size(rgbtmp,4));
            rgbtmp = permute(rgbtmp,[3 2 1]);
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2)*size(rgbtmp,3));
            % ��ʱrgbtmpÿһ�ж�ӦI��һ������������㷨���룬I(m,n)��Ӧ��n+(m-1)*size(I,2)
            I = myAnn24(rgbtmp);
%             I = bpfish24(rgbtmp,r); % ���ﴫ��rֻ��Ϊ����bpfish�����﷽�����Iͼ��Ĵ�С       
        otherwise
    end
%% bpfish��������
%     I = I - double(region(:,:,1));
%     imshow(I)
%% myAnn��������
    I1 = reshape(I,size(I0,2),size(I0,1));
    I1 = I1';
    I1 = I1 - double(region(:,:,1));
    I1(I1<0.5) = 0;
    I1(I1>=0.5) = 1;
    imshow(I1)
end