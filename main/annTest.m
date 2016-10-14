%% Begin ANN Test
clc, close all
neighbors=8;
resize = 1;
[ffile,fpath] = uigetfile({'*.mp4';'*.avi'},'Select the Video file');
fname = [fpath,ffile];
% 读取图层蒙版
region=imread([fname,'_蒙版.png']);
region=imresize(region,resize); %缩放

try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % 读取视频前先判断视频是否已经被读取，可以在调试时减少等待读取视频的时间
    end
catch
    V=VideoReader(fname); 
end

for kk=1:10
    frame = randi(V.NumberOfFrames);
    I0=read(V,frame); % 原始视频帧图像
    I0=imresize(I0,resize); %缩放   
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
            r = double(I(:,:,1)); 
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
%             I = bpfish24(rgbtmp,r); % 这里传入r只是为了在bpfish函数里方便计算I图像的大小       
        otherwise
    end
%% bpfish后续处理
%     I = I - double(region(:,:,1));
%     imshow(I)
%% myAnn后续处理
    I1 = reshape(I,size(I0,2),size(I0,1));
    I1 = I1';
    I1 = I1 - double(region(:,:,1));
    I1(I1<0.5) = 0;
    I1(I1>=0.5) = 1;
    imshow(I1)
end