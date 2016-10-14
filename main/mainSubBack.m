% 基于去背景法的群体运动轨迹追踪算法
% 不使用卡尔曼滤波器，直接用历史平均做预测
% 提高录制视频的帧率能够提高准确度
% 提高视频的背景稳定性有利于减小maxError从而提高在高密度下的正确率
% v1,2014.12.10,能够基本识别
% v2,2014.12.10,修改输出视频帧标示为实际视频对应的帧数,修复有识别到的真实鱼跟踪丢失的情况
% v3,2014.12.10,修复分水岭分割引起的Cshape和fishsegn中的bug
% v4,2014.12.11,修复v3带来的死循环bug,界面小修改,加入maxJump
% v5,2014.12.11,修改信心指数评价方法,在输出txt中增加对应视频中的帧位置
% v6,2014.12.11,增加重新制作视频功能,输出参数中添加了resize,maxJump
% v7,2015.03.16,增加感兴趣区域，初始位置手动选择

% 2015.9.10 需要恰当设置minLight；fishRatio可以不要了
% 2015.9.8  不对最大团进行二次分割。
            % 修改跳跃条件，连续丢帧数越大，可允许的搜索范围越大；（mycost.m）
            % 对于已经丢失的个体也要限定搜索范围，在可允许的搜索范围内的cost一样，之外的Inf （mycost.m）
            % 对于丢失的个体重新出现的，其可能在任何位置是对的，但是如果将成本设置为0则会导致在识别到的个体数量
            % 少于实际个体数量时，丢失的个体占用未丢失个体的位置，而未丢失的个体变成新的丢失个体。
            % 因此，A(ii,:)=9999999999(足够大的数，但不能是Inf)。 （mycost.m）
            % 将亮度小于minLight的点变为全黑，之后再用graythresh函数，效果很好（fishsegn_s1.m）
% 2015.9.7  修改灰度转二值图的level值为0.05（fishsegn_s1.m）
            % 二次分割的图像依然把小于最小面积一半的联通集团去掉 （fishsegn_s2.m）
% 2015.9.6  蚂蚁涂了白色，要识别的个体变成了高亮部分，并且个体大小变得非常小。
            % 如果出现了蚂蚁丢失，设置蚂蚁为静止了。因为很可能是蚂蚁的位置或形态特殊，导致丢失。
            % 正常运动过程中就算有丢失也是断断续续，只有静止了才可能连续多帧丢失。
% 2015.6.11 解决画视频总是置顶的问题
% 2015.4.12 提高对置信度指数的计算准确度
% 2015.4.11 修复背景和前景亮度差异较小时自动灰度转二值图出错的bug，修复了mycost中面积的权重
% 2015.4.10 修复置信度输出一直是Conf0的bug
% 2015.04.09,修改匹配问题为LAP问题，匈牙利算法+自定义成本矩阵，提高易用性。基于mainRgb.m修改。

clc; clear -regexp [^V,^region]; close all; fclose all;
% novideo = 0; % 1为不做结果视频
testmode = 0; % 1为测试模式
testframe = Inf; % 要开始暂停程序的帧
rgb = 0; % 0为利用去背景方法识别
tic
disp('开始识别......')
disp('minLight很重要，切勿忘记恰当设置')
%##########开始参数定义##########
% 还应该根据实际情况调整mycost.m中各个权重大小，目前2,1,2的设置还是很合理的
v95=10; %更保险做法是用v99，95%的情况下鱼速度小于v95,开始可以写大点，如40，测试运行一次，再用testVip.m程序获取真正的v95
minS=50; %用S1,获取方式同v95，一开始可以设置小点，如20，实际只要大于最大误识别的部分就够了
maxS=400; %vmax用S99，获取方式同v95，一开始可以设置大点，如500，设置为略大于真实鱼的最大面积最合适。
           %但识别出的Smax往往是重叠部分的面积，不是真实鱼面积，所以用S95或S99。
% fishRatio = 0.01; % fishRatio代表鱼的面积占整个画面面积的比例，一般略大，但越接近越好。通常该比例小于0.01。
minLight = 10; % 将亮度小于minLight的点变为全黑，之后再用graythresh函数，效果很好（fishsegn_s1.m）
                % 根据实际视频中灰度图里个体核心区域的亮度值来确定minLight，比必须留下的个体像素的最低亮度值略小即可
con=0;
[ffile,fpath] = uigetfile({'*.avi';'*.mp4'},'Select the Video file');
prompt = {'fish number','start frame','frame number','frame interval',...
    'dark fish?','strelSize','resize','maxJump','isabort','novideo'};
name = 'Input Parameters';
numlines = 1;
% minS,maxS应该小点，因为后面做开运算会减小面积
defaultanswer={'20','1','0','1','yes','1','1','1','0','0'}; 
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
novideo = str2double(in{10}); % 1是无视频，其它是有视频
%##########结束参数定义##########
fname=[fpath,ffile];
if testmode==1
    outd=[fpath,'output\test-',ffile,'\'];
else
    outd=[fpath,'output\',ffile,'\'];
end
set_minLight = 0; % 2015.9.10 第一次识别该视频，需要设置minLight
if(~isdir(outd)) 
    mkdir(outd);
    % 2015.9.10 第一次识别该视频，需要设置minLight
    set_minLight = 1;
end

% 读取图层蒙版
% region=imread([fname,'_蒙版.png']);
% region=imresize(region,resize); %缩放



try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % 读取视频前先判断视频是否已经被读取，可以在调试时减少等待读取视频的时间
    end
catch
    V=VideoReader(fname); 
end
try
    region=imread([fname,'_蒙版.png']);
    %region=rgb2gray(region);
catch exception   
 uiwait(msgbox('选择兴趣区域,制作蒙版','Success','modal'));
   
region = roipoly(read(V,100));
region =uint8(~ region.*255);
imwrite(region , [fname,'_蒙版.png']);
end
try
    %如果已经存在背景图片就直接读取，否则就计算背景图片
    back=imread([fname,'_background_',num2str(resize),'.png']);
catch exception
    back=background(V,resize);
    imwrite(back,[fname,'_background_',num2str(resize),'.png']);
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
if V.FrameRate>100
aviobj.FrameRate = 25;
end
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
 frame=sf;
 save('test.mat');
 %用来调节最小面积的，避免出现没有目标的错误
 save('continue.mat','con');
 parameters
 while con==0
     pause(2);
     load continue.mat
 end
 load para.mat
 for frame=sf:interval:sf+(N-1)*interval
    colNow = ceil((frame-sf+1)/interval); %当前帧在处理序列中的序号
    I0=read(V,frame); % 原始视频帧图像
    I0=rgb2gray(I0);
    I0=imresize(I0,resize); %缩放   
    % 灰度模式处理
    if strcmp(darkfish,'yes')
        I = back - I0; %与背景求差；因为背景更亮，所以back-I；反之I-back
    else
        I = I0 - back;
    end
    I=I-region; 
    
%%     % 2015.9.10 设置了minLight，修改了fishsegn_s1后不需要这步了
%     if colNow == 1
%         % 现在I中亮度高的部分是鱼，亮度低(接近0，一般都会小于10)的部分是背景
%         % 但是对于背景和鱼亮度差异小的情况，如果直接自动转换为二值图效果可能会很差 
%         % 所以在传入分割函数前先扩大前景和背景的差异 % 2015.4.11
%         Ivalue = zeros(1,10);
%         Ivaluetmp = reshape(I,1,[]);
%         if fishRatio > 0.5 || fishRatio <= 0
%             disp('fishRatio设置错误，将以默认值0.01进行计算')
%             fishRatio = 0.01;
%         end
%         for tmpii = 1:length(Ivalue)
%             Ivalue(tmpii) = prctile(Ivaluetmp, 100*((1-fishRatio) + (tmpii-1)/length(Ivalue)*fishRatio));
%         end
%         Ivalue = mean(Ivalue);
%     end
%     % 现在的Ivalue基本上可以作为前景和背景亮度的分界线，比Ivalue亮度还小的肯定就是背景了。
%     I(I<Ivalue)=0; % 这可以提高灰度转二值图的正确率
    
    % 2015.9.10 第一次识别该视频，需要设置minLight
    if set_minLight == 1
        imshow(I)
        title('找出前景和背景的分界亮度值，后面将把亮度值低于该阈值的点亮度设为0，以提高后面转二值图的准确率')
        waitInput=1;
        try
            while(waitInput>0)
                pause(str2double(inputdlg('Seconds to wait')));
                waitInput = 0;
            end
        catch
        end
        minLight = str2double(inputdlg('Set minLight'));
        set_minLight = 0;
    end
    
    if colNow<=3
        manual; % 2015-3-17 人工分割确保初始画面无重叠
    end
    
    if testmode==1 && frame>=testframe
        keyboard
    end
   
    
%    imwrite(I,sprintf('fig/%d.png',frame));
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
%                 ax(:,colNow-3),ay(:,colNow-3),len(:,colNow-1),C0(:,1),C0(:,2),S0,v95,maxJump,countLost);
            
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
                % 连续预测了2帧，则这次开始写NaN,预测3帧通常都出错
                Conf (colNow)=Conf(colNow)-20; % 准确度信心指数,每丢失一条鱼减小一定值，连续丢失多帧要减更多信心
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
    fprintf(fid,'   %5.1f       %5.1f    ',frame,Conf(colNow));
    for line=1:num0
        fprintf(fid,'%7.2f  %7.2f  ',px(line,colNow),py(line,colNow));
    end
    fprintf(fid,'\r\n');
    
    if novideo ~= 1 % 无视频模式不画图，不做视频
    if testmode==1 % 测试模式
        imshow(Iseg); hold on
        title('test mode: final segmentation');
    else
            imshow(I0,'border','tight','initialmagnification',100);hold on
    end
    color = hsv(num0);
    for dot=1:num0
        plot(px(dot,colNow), py(dot,colNow),'g.');
        if num0 <= 10
            text(px(dot,colNow), py(dot,colNow),strcat('\color{red}',...
                    num2str(dot)),'FontSize',12);
        end
        linehistory = 20; %画线的长度
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
    'num0','sf','N','interval','darkfish','V',...
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

disp('正常结束')