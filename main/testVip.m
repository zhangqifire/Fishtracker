%% 读视频
clear -regexp [^V | rgbFish | rgbNotFish];
close all; fclose all;clc;
[ffile,fpath] = uigetfile({'*.mp4';'*.avi'},'Select the Video file');
fname = [fpath,ffile];
try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % 读取视频前先判断视频是否已经被读取，可以在调试时减少等待读取视频的时间
    end
catch
    V=VideoReader(fname); 
end
resize=1;
% 读取图层蒙版
region=imread([fname,'_蒙版.png']);
region=imresize(region,resize); %缩放

%% 神经网络识别时，自动定位连续3帧无重叠的位置
num0=1;
strelSize=1;
minS=1;
maxS=300;
neighbors = 8;
tag = 0; % 用于自动定位连续3帧无重叠的位置，即起始识别位置
for ii=100:200
%     frame=randi(V.NumberOfFrames);
    frame = ii;
    I0=read(V,frame); % 原始视频帧图像
    I0=imresize(I0,resize);
    I=I0;
    
    % 这部分是根据红色及邻居颜色梯度信息提取鱼的位置
    % 2015.4.29 加入邻居信息，并将整个图形矩阵输入bpfish.c进行神经网络算法判定
    Itmp = zeros(size(I,1)+4,size(I,2)+4,size(I,3));
    Itmp(3:end-2,3:end-2,:)=double(I);
   
    switch neighbors
        case 0
            rgbtmp = I;
            rgbtmp = permute(rgbtmp,[3 2 1]);
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2)*size(rgbtmp,3));
            % 此时rgbtmp每一列对应I中一个点的神经网络算法输入，I(m,n)对应第n+(m-1)*size(I,2)
            I = myAnn0(rgbtmp);
        case 8
            rgbtmp = zeros(size(I,1),size(I,2),size(I,3),9);
            for jj=0:8
                rgbtmp(:,:,:,jj+1) = Itmp(2+floor(jj/3):end-3+floor(jj/3),2+mod(jj,3):end-3+mod(jj,3),:);
            end 
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2),size(rgbtmp,3)*size(rgbtmp,4));
            rgbtmp = permute(rgbtmp,[3 2 1]);
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2)*size(rgbtmp,3));
            % 此时rgbtmp每一列对应I中一个点的神经网络算法输入，I(m,n)对应第n+(m-1)*size(I,2)
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
            % 此时rgbtmp每一列对应I中一个点的神经网络算法输入，I(m,n)对应第n+(m-1)*size(I,2)
            I = myAnn24(rgbtmp);
        otherwise
    end
    
    I = reshape(I,size(I0,2),size(I0,1)); % 这里reshape方式很关键，不能弄错了。
    I = I';
    I = I - double(region(:,:,1));
    I(I<0.5) = 0;
    I(I>=0.5) = 1;
    I=I-double(region(:,:,1));
    [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s1( I,num0,strelSize,minS,maxS ); %初步识别鱼
    
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

%% 更新蒙版：直接把photoshop处理的蒙版灰色部分改成全黑
Imask=imread([fname,'_蒙版.png']);
for ii=1:size(Imask,1)
    for jj=1:size(Imask,2)
        if Imask(ii,jj,1)>0 || Imask(ii,jj,2)>0 || Imask(ii,jj,3)>0
            Imask(ii,jj,:)=[255,255,255];
        end
    end
end
imwrite(Imask,[fname,'_蒙版.png'],'png');
imshow(Imask)

%% 鱼外形的基本统计信息
% 对于新的视频，先找一段好识别的进行识别，利用这些数据进行初步统计
% clear -regexp [^V]
% load '40pie100-400.avi_sf=1_ef=301.mat';
% load 'VID_20130707_094142_草金鱼旋转1.mp4_sf=1_N=300.mat';
vmax=max(max(sqrt(vx.^2+vy.^2))); % 46.1
v1=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),1); % 分位点数据 0 
v5=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),5); % 分位点数据 0
v99=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),99); % 分位点数据 39.6 
v95=prctile(reshape(sqrt(vx.^2+vy.^2),1,[]),95); % 分位点数据 29.7
Smax=max(max(S)); % 212
S1=prctile(reshape(S,1,[]),1); % 分位点数据 24.5
S5=prctile(reshape(S,1,[]),5); % 分位点数据 28
S99=prctile(reshape(S,1,[]),99); % 分位点数据 176
S95=prctile(reshape(S,1,[]),95); % 分位点数据 170.5
% Ecmin=min(min(Ec));
% Ecmax=max(max(Ec));
% Ec1=prctile(reshape(Ec,1,[]),1); % 分位点数据
% Ec99=prctile(reshape(Ec,1,[]),99); % 分位点数据
% Ec5=prctile(reshape(Ec,1,[]),5); % 分位点数据
% Ec95=prctile(reshape(Ec,1,[]),95); % 分位点数据
% lenmin=min(min(len));
% lenmax=max(max(len));
% len1=prctile(reshape(len,1,[]),1); % 分位点数据
% len99=prctile(reshape(len,1,[]),99); % 分位点数据
% widmin=min(min(wid));
% widmax=max(max(wid));
% wid1=prctile(reshape(wid,1,[]),1); % 分位点数据
% wid99=prctile(reshape(wid,1,[]),99); % 分位点数据
