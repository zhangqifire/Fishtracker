%% 识别过程中为了方便没有做视频，通过本程序可以利用识别数据制作识别结果视频
% 读取参数
clc; clear all; close all; fclose all;
tic
[ffile,fpath] = uigetfile({'*.avi';'*.mp4'},'Select the Video file');
fname=[fpath,ffile];
V=VideoReader(fname);
[fdata,fpath] = uigetfile([fpath,'*.mat'],'Select the Data file'); %保持上一步选择的路径信息
fdata=[fpath,fdata];
load(fdata);
[fdata,fpath] = uigetfile([fpath,'*.mat'],'Select the Parameters file');
fdata=[fpath,fdata];
load(fdata);
if ~exist('resize','var')
    resize = 1;
end
%% 开始重画视频
fname=[fpath,ffile]; %文件名去后缀
% 初始化输出视频
aviobj = VideoWriter(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N)),'MPEG-4');
aviobj.FrameRate = 25;
open(aviobj);

for frame=sf+2*interval:interval:sf+(N-1)*interval 
    colNow = ceil((frame-sf+1)/interval); %当前帧在处理序列中的序号
    I0=read(V,frame); % 原始视频帧图像
    I0=imresize(I0,resize); %缩放   
    C0(:,1)=px(:,colNow); %t1,t2时刻x坐标，图像中的x,y轴
    C0(:,2)=py(:,colNow); %t1,t2时刻y坐标
    
    imshow(I0,'border','tight','initialmagnification',100);
%     imshow(I0);
    hold on
    color = hsv(num0);
    for dot=1:num0
        plot(px(dot,colNow), py(dot,colNow),'g.');
        if num0 <= 10
            text(px(dot,colNow), py(dot,colNow),strcat('\color{red}',...
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

%%
close(aviobj);
close('all');
fclose('all');
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.mp4'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'.mp4'))
tElapsed = toc;
disp(['重新画图总时间：',num2str(tElapsed),'s 平均每帧：',num2str(tElapsed/N),'s'])
