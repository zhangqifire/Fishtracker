% 基于彩色信息的鱼群分割及轨迹追踪算法
% 针对'VID_20130707_094142_草金鱼旋转1.mp4’这类视频背景不稳定，灰度信息不稳定，
% 现有识别软件及方法效果差的情况开发。
% step 1. 分割：先将所有画面的鱼分割出来，鱼看上去都是红色的，肉眼很容易区分，
% 根据RGB值判定任意色块是否属于鱼。使用规则判定和人工神经网络算法判定。
% setp 2. 连接：自定义成本函数（与速度、加速度、面积和历史信息相关），用匈牙利算法匹配。

% 在使用彩色信息识别时，grayThreshold等参数用不上，用灰度识别时用得上，改为程序自动计算；
% 鱼的大小不一样，所以原来的minS,maxS需要由比例改为实际像素个数，直接根据前几帧数据得出；
% 现在的预测算法确定只需要3帧历史信息，history参数也不需要了

% 2015.9.9  选择有效数据的外框，减少计算负担
% 2015.9.8  不对最大团进行二次分割。
            % 修改跳跃条件，连续丢帧数越大，可允许的搜索范围越大；（mycost.m）
            % 对于已经丢失的个体也要限定搜索范围，在可允许的搜索范围内的cost一样，之外的Inf （mycost.m）
            % 对于丢失的个体重新出现的，其可能在任何位置是对的，但是如果将成本设置为0则会导致在识别到的个体数量
            % 少于实际个体数量时，丢失的个体占用未丢失个体的位置，而未丢失的个体变成新的丢失个体。
            % 因此，A(ii,:)=9999999999(足够大的数，但不能是Inf)。 （mycost.m）
            % 将亮度小于minLight的点变为全黑，之后再用graythresh函数，效果很好（fishsegn_s1.m）
% 2015.9.7  修改灰度转二值图的level值为0.05（fishsegn_s1.m）
            % 二次分割的图像依然把小于最小面积一半的联通集团去掉 （fishsegn_s2.m）
% 2015.6.11 解决画视频总是置顶的问题
% 2015.5.1  调用matlab自带神经网络工具包训练结果myAnn0/8/24.m进行判定
% 2015.4.29 调用新的bpfish.c进行判定
% 2015.4.12 提高对置信度指数的计算准确度
% 2015.4.10 修复置信度输出一直是Conf0的bug;
%          改为完全用神经网络确定，因为加入规则的话每次换视频除了要重新训练神经网络还要改判定规则，很麻烦
%          并且神经网络的分界线确定为(y<=0.5)?0.0:1.0，这个可以训练做到
% 2015/4/9 修复mycost程序NaN数据导致的程序崩溃，增加不做视频模式，修改最大连续预测帧数为2
%          修复信心指数的一个错误
% 2015/4/8 修改连接算法为匈牙利算法匹配连接；修改参数设置；修改丢帧处理
% 2015/4/2 先只修改识别程序，看看颜色分割是否可行
clc; clear -regexp [^V]; close all; fclose all;
% novideo = 0; % 1为不做结果视频
testmode = 0; % 1为测试模式
testframe = Inf; % 要开始暂停程序的帧
rgb = 1; % 1为利用rgb信息识别
tic
disp('远程桌面操作时不可以最小化或断开连接，否则制作的视频会黑屏')
%##########开始参数定义##########
% 还应该根据实际情况调整mycost.m中各个权重大小，目前2,1,2的设置还是很合理的
v95=10; %95%的情况下鱼速度小于v95,开始可以写大点，如40，测试运行一次，再用testVip.m程序获取真正的v95
minS=200; %获取方式同v95，一开始可以设置小点，如20
maxS=500; %获取方式同v95，一开始可以设置大点，如500
minLight = 0.5; % 2015.9.9 用神经网络识别的特殊处理，因为传入fishsegn_s1.m的就已经是二值图了，
                % 所以不需要转二值图，直接处理联通集团即可。
                % 将亮度小于minLight的点变为全黑，之后再用graythresh函数，效果很好（fishsegn_s1.m）
                % 根据实际视频中灰度图里个体核心区域的亮度值来确定minLight，比必须留下的个体像素的最低亮度值略小即可

[ffile,fpath] = uigetfile({'*.mp4';'*.avi'},'Select the Video file');
prompt = {'fish number','start frame','frame number','frame interval',...
    'dark fish?','strelSize','resize','maxJump','isabort','isbp','novideo','neighbors'};
name = 'Input Parameters';
numlines = 1;
% minS,maxS应该小点，因为后面做开运算会减小面积
defaultanswer={'1','27000','0','1','yes','0','1','1','0','1','0','8'}; 
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
in = inputdlg(prompt,name,numlines,defaultanswer,options);
num0 = str2double(in{1}); % 已知鱼数量
sf = str2double(in{2}); % 开始识别的帧，鱼在此帧及下interval帧中不能有重叠
N = str2double(in{3});  % 识别的帧数,如果为0则识别整个视频
interval = str2double(in{4}); % 数据处理帧间隔
darkfish = in{5}; % 鱼比背景更暗为yes，否则为no
strelSize = str2double(in{6}); % 形态学操作结构元素大小，一般可以取4
resize = str2double(in{7}); % 缩放比例
maxJump = str2double(in{8}); % 最大跳跃系数，基本上相当于误差范围
isabort = str2double(in{9}); % 为1表示如果同时出现了2条鱼为NaN则结束程序
isbp = str2double(in{10}); % 为1表示采用神经网络识别
novideo = str2double(in{11}); % 1是无视频，其它是有视频
neighbors = str2double(in{12}); % 0,8,24个邻居
%##########结束参数定义##########
fname=[fpath,ffile];
if testmode==1
    outd=[fpath,'output\test-',ffile,'\'];
else
    outd=[fpath,'output\',ffile,'\'];
end
if(~isdir(outd)) 
    mkdir(outd); 
end

% 读取图层蒙版
region=imread([fname,'_蒙版.png']);


% 2015.9.9 选择有效数据的外框，减少计算负担
imshow(region)
title('现在选择有效数据的外框','Color','b')
[~, rect] = imcrop(); % 交互式选择,[xmin ymin width height] 
rmin=round(rect(2));rmax=round(rmin+rect(4));cmin=round(rect(1));cmax=round(cmin+rect(3)); %外
close gcf

try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % 读取视频前先判断视频是否已经被读取，可以在调试时减少等待读取视频的时间
    end
catch
    V=VideoReader(fname); 
end

if N==0 || sf+(N-1)*interval > V.NumberOfFrames
    % 如果N为0或者要识别的最大帧超出了视频总帧数，则N设置为最大值
    N = (V.NumberOfFrames - sf) / interval + 1;
end
% 初始化输出坐标文件
fname=[outd,ffile]; %文件名去后缀
fid=strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.txt');
fid2=strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_lost.txt'); %丢帧信息
fid2=fopen(fid2,'w'); 

% 写入输出文件
fid=fopen(fid,'w');
% % 写入标题' Confidence   px1    py1   px2   py2...'
fprintf(fid, ' FramePos  Confidence ');
for line=1:num0
    fprintf(fid,['   px',int2str(line),'      py',int2str(line),'   ']);
end
fprintf(fid,'\r\n');

if novideo ~= 1 % 无视频模式不画图，不做视频
% 初始化输出视频
aviobj = VideoWriter(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N)),'MPEG-4');
aviobj.FrameRate = V.FrameRate;
aviobj.Quality = 100;
open(aviobj);  
end

% 初始化要记录的变量
px=NaN(num0,N); % x坐标位置
py=px; % y坐标位置
vx=px; % x方向速度
vy=px; % y方向速度
ax=px; % x方向加速度
ay=py; % y方向加速度
len=px; % 椭圆拟合得到的体长
wid=px; % 椭圆拟合得到的体宽
Ec=px; % 椭圆离心率
S=px; % 鱼的面积
Conf=NaN(1,N); % 准确度信心指数

countLost = zeros(num0,1); % 初始化丢帧计数器
abort = 0; % 初始化终止程序标识符
target_last = zeros(num0,1); % 记录上一步的匹配结果
% 冷启动，刚开始处理的帧
% 冷启动时要求鱼无重叠，需要三帧画面，可以用手工分割
for frame=sf:interval:sf+(N-1)*interval
    colNow = ceil((frame-sf+1)/interval); %当前帧在处理序列中的序号
    I0=read(V,frame); % 原始视频帧图像
    I0=imresize(I0,resize); %缩放   
    I=I0;
%     I=I0-region;
    
    if colNow<=3
        manual; % 2015-3-17 人工分割确保初始画面无重叠
    end
    
    % 这部分是根据红色及邻居颜色梯度信息提取鱼的位置
    % 2015.4.29 加入邻居信息，并将整个图形矩阵输入bpfish.c进行神经网络算法判定
    % 2015.5.1  改为用matlab自带ANN工具箱训练，因为速度更快
    % 2015.9.9
    I = I(rmin:rmax,cmin:cmax,:); % 有效数据外框，只需要处理框内的数据即可
    Idata = I;
    
    Itmp = zeros(size(I,1)+4,size(I,2)+4,size(I,3));
    Itmp(3:end-2,3:end-2,:)=double(I);
    switch neighbors
        case 0
            rgbtmp = I;
            rgbtmp = permute(rgbtmp,[3 2 1]);
            rgbtmp = reshape(rgbtmp,size(rgbtmp,1),size(rgbtmp,2)*size(rgbtmp,3));
            % 此时rgbtmp每一列对应I中一个点的神经网络算法输入，I(m,n)对应第n+(m-1)*size(I,2)
            rgbtmp = double(rgbtmp);
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
    
    % 2015.9.9
    I = reshape(I,size(Idata,2),size(Idata,1)); % 这里reshape方式很关键，不能弄错了。
%     I = reshape(I,size(I0,2),size(I0,1)); % 这里reshape方式很关键，不能弄错了。
    I = I';
    % 2015.9.9
    Idata = zeros(size(I0,1),size(I0,2)); % 把图全变为0，再把有效数据部分替换为识别结果
    Idata(rmin:rmax,cmin:cmax) = I; % 有效数据外框，只需要处理框内的数据即可
    I = Idata;
%     I=I-double(region(:,:,1));
%     I(I<0.5) = 0;
%     I(I>=0.5) = 1;
    I = im2bw(I) - im2bw(region);
    [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s1( I,num0,strelSize,minS,maxS,minLight ); %初步识别鱼
    if colNow > 3
        % 已经过了初始化的阶段，则考虑根据匹配结果二次分割
        % 先预匹配一次，根据结果对范围内大于面积均值的集团进行分割
        % 求成本矩阵 source行，target列         
        A = mycost(px(:,colNow-1),py(:,colNow-1),vx(:,colNow-2),vy(:,colNow-2),S(:,colNow-1), ...
            ax(:,colNow-3),ay(:,colNow-3),len(:,colNow-1),C0(:,1),C0(:,2),S0,v95,maxJump,countLost);
        % 将一条鱼识别为多条鱼没什么问题，但是识别不出鱼会比较难办，所以要尽量识别出鱼来。
        % 同一个杂质在连续多帧中被误识别出的概率较小，所以可允许一定的误识别。
        target_indices = assignmentoptimal(A); % 此处A为成本，越小越好
        ind = find(target_indices==0); % 查找没匹配上的
        while ~isempty(ind)
            % 如果存在没有匹配上的，则对处于该个体可能范围内的最大个体进行分割            
            indtmp = ind(1); % 处理indtmp号鱼
            ind(1) = []; % 将ind(1)移除
            % 上一步的位置为圆心，以v95+maxJump*len为半径
            pxtmp = px(indtmp,colNow-1); %x坐标，图像中的x,y轴
            pytmp = py(indtmp,colNow-1); %y坐标
            distmp = v95 + maxJump * len(indtmp,colNow-1); % 可能范围大小
            disC0 = sqrt((pxtmp-C0(:,1)).^2 + (pytmp-C0(:,2)).^2); % 与C0各点距离
            C0indtmp = find(disC0<=distmp); % 需要分割的连通集团序号
            if ~isempty(C0indtmp)
                % 必须有需要分割的团才分割，不然要出错
                [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s2( Iseg,C0indtmp,Conf0,num0,minS );
            end
        end
    end
    Conf(colNow)=Conf0; % 准确度信心指数,与所有鱼相关
    
    if colNow==1
        px(:,colNow)=C0(:,1); %x坐标，图像中的x,y轴
        py(:,colNow)=C0(:,2); %y坐标
        len(:,colNow)=len0; % 椭圆拟合得到的体长
        wid(:,colNow)=wid0; % 椭圆拟合得到的体宽
        Ec(:,colNow)=Ec0; % 椭圆离心率
        S(:,colNow)=S0; % 鱼的面积
        Conf(colNow)=Conf0; % 准确度信心指数
    else
        if colNow<=3
            source = [px(:,colNow-1),py(:,colNow-1)];
            target = [C0(:,1),C0(:,2)];
            % 如果是刚开始，就用hungarianlinker.m进行匹配
            % 开始几步不设置限制条件的最优匹配，所以要求识别准确
            target_indices = hungarianlinker(source, target); 
        else
            % 求成本矩阵 source行，target列         
            A = mycost(px(:,colNow-1),py(:,colNow-1),vx(:,colNow-2),vy(:,colNow-2),S(:,colNow-1), ...
                ax(:,colNow-3),ay(:,colNow-3),len(:,colNow-1),C0(:,1),C0(:,2),S0,v95,maxJump,countLost);
              %%%%%%% test for A %%%%%%%%%%%%
%             frame
%             [A,Av,Aa,As] = mycost(px(:,colNow-1),py(:,colNow-1),vx(:,colNow-2),vy(:,colNow-2),S(:,colNow-1), ...
%                ax(:,colNow-3),ay(:,colNow-3),len(:,colNow-1),C0(:,1),C0(:,2),S0,v95,maxJump,countLost);
            
            % 将一条鱼识别为多条鱼没什么问题，但是识别不出鱼会比较难办，所以要尽量识别出鱼来。
            % 同一个杂质在连续多帧中被误识别出的概率较小，所以可允许一定的误识别。
            target_indices = assignmentoptimal(A); % 此处A为成本，越小越好
            
            %%%%% 根据成本矩阵更新置信度 %%%%%%%%%%%
            % 如果这次匹配的结果与上次不一样，且成本矩阵存在相近的情况，则降低置信度
            if ~isequal(target_indices, target_last)
                Aminr = min(A,[],2); % A中每一行的最小值
                B = ( A - repmat(Aminr,1,size(A,2))) ./ repmat(Aminr,1,size(A,2)) ; % B是成本函数的相似度
                if sum(sum(B<0.1,2) >= 2) >= 2
                    Conf(colNow)=Conf(colNow)-30; % 如果有2个个体有2个以上的可能匹配，则很容易出错
                end
            end
            
            if testmode==1 && frame>=testframe
                disp(A)
                disp(target_indices')
                keyboard
            end
            target_last = target_indices; % 更新上一步的匹配结果
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % 利用帧间配对信息更新识别数据        
        indyes = find(target_indices~=0); % 查找匹配上的
        while ~isempty(indyes)
            % 对于没有丢帧的鱼，则正常更新状态
            indtmp = indyes(1); % 处理indtmp号鱼
            indyes(1) = []; % 将indyes(1)移除
            countLost(indtmp) = 0; % 重置丢帧计数器
            px(indtmp,colNow)=C0(target_indices(indtmp),1); %x坐标，图像中的x,y轴
            py(indtmp,colNow)=C0(target_indices(indtmp),2); %y坐标
            len(indtmp,colNow)=len0(target_indices(indtmp)); % 椭圆拟合得到的体长
            wid(indtmp,colNow)=wid0(target_indices(indtmp)); % 椭圆拟合得到的体宽
            Ec(indtmp,colNow)=Ec0(target_indices(indtmp)); % 椭圆离心率
            S(indtmp,colNow)=S0(target_indices(indtmp)); % 鱼的面积
            vx(indtmp,colNow-1) = px(indtmp,colNow) - px(indtmp,colNow-1); % 速度
            vy(indtmp,colNow-1) = py(indtmp,colNow) - py(indtmp,colNow-1);
            if colNow>2
                ax(indtmp,colNow-2) = px(indtmp,colNow) + px(indtmp,colNow-2) - 2 * px(indtmp,colNow-1); % 加速度
                ay(indtmp,colNow-2) = py(indtmp,colNow) + py(indtmp,colNow-2) - 2 * py(indtmp,colNow-1); 
            end
        end
       
        ind = find(target_indices==0); % 查找没匹配上的
        while ~isempty(ind)
            % 如果出现鱼丢失，将丢失鱼的信息用历史信息预测填充
            disp(['出现了丢帧，帧位置为：',num2str(frame)])
            fprintf(fid2,['出现了丢帧，帧位置为：',num2str(frame),'\r\n']);
            % 加入自适应的速度大小限制
            vel = sqrt((vx(:,colNow - 2).^2 + vy(:,colNow -2).^2)); % t-2时刻速度大小
            vxtmp = vx(:,colNow-2)*2 - vx(:,colNow-3);
            vytmp = vy(:,colNow-2)*2 - vy(:,colNow-3);
            vtmp = sqrt(vxtmp.^2 + vytmp.^2); % t-1时刻速度大小
            pxtmp = px(:,colNow-1) + vxtmp .* vel ./ vtmp;
            pytmp = py(:,colNow-1) + vytmp .* vel ./ vtmp;
            
            indtmp = ind(1); % 处理indtmp号鱼
            ind(1) = []; % 将ind(1)移除
            countLost(indtmp) = countLost(indtmp) + 1;
            if countLost(indtmp) > 2
                % 连续预测了若干帧，则这次开始写NaN,预测3帧通常都出错
                Conf(colNow)=Conf(colNow)-20; % 准确度信心指数,每丢失一条鱼减小一定值，连续丢失多帧要减更多信心
                px(indtmp,colNow)=NaN; %x坐标，图像中的x,y轴
                py(indtmp,colNow)=NaN; %y坐标
                len(indtmp,colNow)=NaN; % 椭圆拟合得到的体长
                wid(indtmp,colNow)=NaN; % 椭圆拟合得到的体宽
                Ec(indtmp,colNow)=NaN; % 椭圆离心率
                S(indtmp,colNow)=NaN; % 鱼的面积
                vx(indtmp,colNow-1) = NaN; % 速度
                vy(indtmp,colNow-1) = NaN;
                if colNow>2
                    ax(indtmp,colNow-2) = NaN; % 加速度
                    ay(indtmp,colNow-2) = NaN; 
                end
            else
                Conf(colNow)=Conf(colNow)-10; % 准确度信心指数,每丢失一条鱼减小一定值
                px(indtmp,colNow)=pxtmp(indtmp); %x坐标，图像中的x,y轴
                py(indtmp,colNow)=pytmp(indtmp); %y坐标
                len(indtmp,colNow)=len(indtmp,colNow-1); % 椭圆拟合得到的体长
                wid(indtmp,colNow)=wid(indtmp,colNow-1); % 椭圆拟合得到的体宽
                Ec(indtmp,colNow)=Ec(indtmp,colNow-1); % 椭圆离心率
                S(indtmp,colNow)=S(indtmp,colNow-1); % 鱼的面积
                vx(indtmp,colNow-1) = px(indtmp,colNow) - px(indtmp,colNow-1); % 速度
                vy(indtmp,colNow-1) = py(indtmp,colNow) - py(indtmp,colNow-1);
                if colNow>2
                    ax(indtmp,colNow-2) = px(indtmp,colNow) + px(indtmp,colNow-2) - 2 * px(indtmp,colNow-1); % 加速度
                    ay(indtmp,colNow-2) = py(indtmp,colNow) + py(indtmp,colNow-2) - 2 * py(indtmp,colNow-1); 
                end
            end
        end

        if sum(isnan(px(:,colNow)))>=2 && isabort==1
            % 如果有2个或更多个体已经彻底丢失则终止程序
            abort = 1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    end
    
    % % 写入识别结果
    fprintf(fid,'   %5.1f       %5.1f    ',frame,Conf(colNow)); % 2015.4.10 修复置信度输出一直是Conf0的bug
    for line=1:num0
        fprintf(fid,'%7.2f  %7.2f  ',px(line,colNow),py(line,colNow));
    end
    fprintf(fid,'\r\n');
    
    if novideo ~= 1 % 无视频模式不画图，不做视频
    if testmode==1 % 测试模式
        imshow(Iseg); hold on
        title('test mode: final segmentation');
    else
        imshow(I0);hold on
    end
    color = hsv(num0);
    for dot=1:num0
        plot(px(dot,colNow), py(dot,colNow),'g.');
        if num0 <= 10
            text(px(dot,colNow), py(dot,colNow),strcat('\leftarrow \color{red}',...
                    num2str(dot)),'FontSize',12);
        end
        linehistory = 30; %画线的长度
        if colNow<=linehistory
            plot(px(dot,1:colNow),py(dot,1:colNow),'Color',color(dot,:))
        elseif colNow>linehistory
            plot(px(dot,colNow-linehistory:colNow),py(dot,colNow-linehistory:colNow), ...
                'Color',color(dot,:))
        end
    end
    text(3,20,['N = ',int2str(size(C0,1)),'  Frame = ',int2str(frame)],...
        'fontsize',10,'color',[1 1 1],'BackgroundColor',[.5 .5 .5]);
    hold off
    drawnow;      
%     frameout = getframe;
    % 2015.6.11 解决画视频总是置顶的问题
    frameout = im2frame(zbuffer_cdata(gcf));
    writeVideo(aviobj,frameout);  
    end
    if abort==1
        break
    end
end
%% 保存最终结果到MAT文件
px(:,colNow+1:end)=[]; % 清理未识别部分，减小文件大小
py(:,colNow+1:end)=[];
ax(:,colNow+1:end)=[];
ay(:,colNow+1:end)=[];
vx(:,colNow+1:end)=[];
vy(:,colNow+1:end)=[];
len(:,colNow+1:end)=[];
wid(:,colNow+1:end)=[];
Ec(:,colNow+1:end)=[];
S(:,colNow+1:end)=[];
Conf(colNow+1:end)=[];
save(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.mat'), ...
    'px','py','vx','vy','ax','ay','len','wid','Ec','S','Conf');

tElapsed = toc % 显示并保存程序运行所花时间
fprintf(fid2,['tElapsed = ',num2str(tElapsed),' s\r\n']);
% 保存参数到MAT文件
save(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_par.mat'), ...
    'num0','sf','N','interval','darkfish','V','isbp',...
    'strelSize','tElapsed','resize','maxJump','isabort','v95','minS','maxS','minLight');
if novideo ~= 1 % 无视频模式不画图，不做视频
close(aviobj);
end
close('all');
fclose('all');
se = frame; %实际识别的最后一帧
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.txt'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'.txt'))
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.mat'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'.mat'))
if novideo ~= 1 % 无视频模式不画图，不做视频
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.mp4'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'.mp4'))
end
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_lost.txt'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'_lost.txt'))
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_par.mat'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'_par.mat'))
load handel.mat; %程序结束，播放提示音
sound(y,2*Fs);