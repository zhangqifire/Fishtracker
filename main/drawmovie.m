%% ʶ�������Ϊ�˷���û������Ƶ��ͨ���������������ʶ����������ʶ������Ƶ
% ��ȡ����
clc; clear all; close all; fclose all;
tic
[ffile,fpath] = uigetfile({'*.avi';'*.mp4'},'Select the Video file');
fname=[fpath,ffile];
V=VideoReader(fname);
[fdata,fpath] = uigetfile([fpath,'*.mat'],'Select the Data file'); %������һ��ѡ���·����Ϣ
fdata=[fpath,fdata];
load(fdata);
[fdata,fpath] = uigetfile([fpath,'*.mat'],'Select the Parameters file');
fdata=[fpath,fdata];
load(fdata);
if ~exist('resize','var')
    resize = 1;
end
%% ��ʼ�ػ���Ƶ
fname=[fpath,ffile]; %�ļ���ȥ��׺
% ��ʼ�������Ƶ
aviobj = VideoWriter(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N)),'MPEG-4');
aviobj.FrameRate = 25;
open(aviobj);

for frame=sf+2*interval:interval:sf+(N-1)*interval 
    colNow = ceil((frame-sf+1)/interval); %��ǰ֡�ڴ��������е����
    I0=read(V,frame); % ԭʼ��Ƶ֡ͼ��
    I0=imresize(I0,resize); %����   
    C0(:,1)=px(:,colNow); %t1,t2ʱ��x���꣬ͼ���е�x,y��
    C0(:,2)=py(:,colNow); %t1,t2ʱ��y����
    
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
        linehistory = 30; %���ߵĳ���
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
    % 2015.6.11 �������Ƶ�����ö�������
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
disp(['���»�ͼ��ʱ�䣺',num2str(tElapsed),'s ƽ��ÿ֡��',num2str(tElapsed/N),'s'])
