%% ����Ƶ
clear -regexp [^V | rgbFish | rgbNotFish];
close all; fclose all;clc;
[ffile,fpath] = uigetfile({'*.mp4';'*.avi'},'Select the Video file');
fname = [fpath,ffile];
try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % ��ȡ��Ƶǰ���ж���Ƶ�Ƿ��Ѿ�����ȡ�������ڵ���ʱ���ٵȴ���ȡ��Ƶ��ʱ��
    end
catch
    V=VideoReader(fname); 
end
resize=1;
% ��ȡͼ���ɰ�
region=imread([fname,'_�ɰ�.png']);
region=imresize(region,resize); %����

%% ������ʶ��ʱ���Զ���λ����3֡���ص���λ��
num0=1;
strelSize=1;
minS=1;
maxS=300;
neighbors = 8;
tag = 0; % �����Զ���λ����3֡���ص���λ�ã�����ʼʶ��λ��
for ii=100:200
%     frame=randi(V.NumberOfFrames);
    frame = ii;
    I0=read(V,frame); % ԭʼ��Ƶ֡ͼ��
    I0=imresize(I0,resize);
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
        otherwise
    end
    
    I = reshape(I,size(I0,2),size(I0,1)); % ����reshape��ʽ�ܹؼ�������Ū���ˡ�
    I = I';
    I = I - double(region(:,:,1));
    I(I<0.5) = 0;
    I(I>=0.5) = 1;
    I=I-double(region(:,:,1));
    [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s1( I,num0,strelSize,minS,maxS ); %����ʶ����
    
%     figure;
    imshow(Iseg)
    title (['num=',num2str(length(S0))])
    drawnow
    if length(S0)==num0
        tag=tag+1;
        if tag==3
            break
        end
    else
        tag=0;
    end
end
disp(ii-2)

%% �����ɰ棺ֱ�Ӱ�photoshop������ɰ��ɫ���ָĳ�ȫ��
Imask=imread([fname,'_�ɰ�.png']);
for ii=1:size(Imask,1)
    for jj=1:size(Imask,2)
        if Imask(ii,jj,1)>0 || Imask(ii,jj,2)>0 || Imask(ii,jj,3)>0
            Imask(ii,jj,:)=[255,255,255];
        end
    end
end
imwrite(Imask,[fname,'_�ɰ�.png'],'png');
imshow(Imask)

%% �����εĻ���ͳ����Ϣ
% �����µ���Ƶ������һ�κ�ʶ��Ľ���ʶ��������Щ���ݽ��г���ͳ��
% clear -regexp [^V]
% load '40pie100-400.avi_sf=1_ef=301.mat';
% load 'VID_20130707_094142_�ݽ�����ת1.mp4_sf=1_N=300.mat';
vmax=max(max(sqrt(vx.^2+vy.^2))); % 46.1
v1=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),1); % ��λ������ 0 
v5=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),5); % ��λ������ 0
v99=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),99); % ��λ������ 39.6 
v95=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),95); % ��λ������ 29.7
Smax=max(max(S)); % 212
S1=prctile(reshape(S,1,[]),1); % ��λ������ 24.5
S5=prctile(reshape(S,1,[]),5); % ��λ������ 28
S99=prctile(reshape(S,1,[]),99); % ��λ������ 176
S95=prctile(reshape(S,1,[]),95); % ��λ������ 170.5
% Ecmin=min(min(Ec));
% Ecmax=max(max(Ec));
% Ec1=prctile(reshape(Ec,1,[]),1); % ��λ������
% Ec99=prctile(reshape(Ec,1,[]),99); % ��λ������
% Ec5=prctile(reshape(Ec,1,[]),5); % ��λ������
% Ec95=prctile(reshape(Ec,1,[]),95); % ��λ������
% lenmin=min(min(len));
% lenmax=max(max(len));
% len1=prctile(reshape(len,1,[]),1); % ��λ������
% len99=prctile(reshape(len,1,[]),99); % ��λ������
% widmin=min(min(wid));
% widmax=max(max(wid));
% wid1=prctile(reshape(wid,1,[]),1); % ��λ������
% wid99=prctile(reshape(wid,1,[]),99); % ��λ������
