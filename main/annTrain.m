%% initial
%clc; clear -regexp [^V];
global neighbors I;
neighbors = 8; % 0, 8, 24

resize = 1;
 %rgbFish = [];
 %rgbNotFish = [];
%load 'rgbFish.mat';
%load 'rgbNotFish.mat';

%% open video
[ffile,fpath] = uigetfile({'*.mp4';'*.avi'},'Select the Video file');
fname = [fpath,ffile];
outd = [fpath,'annGUI\',ffile,'\'];
if(~isdir(outd)) 
    mkdir(outd); 
end
try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % ��ȡ��Ƶǰ���ж���Ƶ�Ƿ��Ѿ�����ȡ�������ڵ���ʱ���ٵȴ���ȡ��Ƶ��ʱ��
    end
catch
    V=VideoReader(fname); 
end

% ��ȡͼ���ɰ�
region=imread([fname,'_�ɰ�.png']);
region=imresize(region,resize); %����

%% choose fish data
go = 'yes';
while(strcmp(go,'yes')==1)
    frame = randi(V.NumberOfFrames);
    I = read(V,frame); % ԭʼ��Ƶ֡ͼ��
    I = imresize(I,resize); %����   
    %I = I-region;
    imshow(I);
    title('choose some fish data points and export cursor\_info')
    cursor_info=1;
    try
        while(cursor_info==1)
            pause(1);
        end
    catch
    end
    
    % ��ȡ���ݵ����꼰rgb��Ϣ
    pos=cat(1,cursor_info.Position); 
        
    % ����pos��ȡ��ɫѵ����
    posN = size(pos,1);
    switch neighbors
        case 0
            rgbtmp = zeros(posN,3);
        case 8
            rgbtmp = zeros(posN,3*9);
        case 24
            rgbtmp = zeros(posN,3*25);   
        otherwise
    end
    for ii=1:posN
        x = pos(ii,2); % pos��ͼ���е�x,y;��Ӧ�������x����y
        y = pos(ii,1);
        % ��ȡ���RGB��Ϣ
        switch neighbors
            case 0
                rgbtmp(ii,:) = I(x,y,:); 
            case 8
                for jj=0:8
                    rgbtmp(ii,jj*3+1:jj*3+3) = I(x-1+floor(jj/3),y-1+mod(jj,3),:);
                end
            case 24
                for jj=0:24
                    rgbtmp(ii,jj*3+1:jj*3+3) = I(x-2+floor(jj/5),y-2+mod(jj,5),:);
                end                
            otherwise
        end
    end
    
    % �ϲ�rgbtmp������ɫѵ�����ݼ�
    rgbFish=[rgbFish;rgbtmp];
    
    prompt1 = {'go on?'};
    name1 = 'Segmenting manual';
    numlines1 = 1;
    defaultanswer1={'yes'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    in = inputdlg(prompt1,name1,numlines1,defaultanswer1,options);
    try
        go = in{1};
    catch
        go = 'no';
    end
end
save([outd,'rgbFish.mat'],'rgbFish');

%% choose background (not fish) data
go = 'yes';
while(strcmp(go,'yes')==1)
    frame = randi(V.NumberOfFrames);
    I = read(V,frame); % ԭʼ��Ƶ֡ͼ��
    I = imresize(I,resize); %����   
%    I = I-region;
    imshow(I);
    title('choose some background data points and export cursor\_info')
    cursor_info=1;
    try
        while(cursor_info==1)
            pause(1);
        end
    catch
    end
    
    % ��ȡ���ݵ����꼰rgb��Ϣ
    pos=cat(1,cursor_info.Position); 
        
    % ����pos��ȡ��ɫѵ����
    posN = size(pos,1);
    switch neighbors
        case 0
            rgbtmp = zeros(posN,3);
        case 8
            rgbtmp = zeros(posN,3*9);
        case 24
            rgbtmp = zeros(posN,3*25);   
        otherwise
    end
    for ii=1:posN
        x = pos(ii,2); % pos��ͼ���е�x,y;��Ӧ�������x����y
        y = pos(ii,1);
        % ��ȡ���RGB��Ϣ
        switch neighbors
            case 0
                rgbtmp(ii,:) = I(x,y,:); 
            case 8
                for jj=0:8
                    rgbtmp(ii,jj*3+1:jj*3+3) = I(x-1+floor(jj/3),y-1+mod(jj,3),:);
                end
            case 24
                for jj=0:24
                    rgbtmp(ii,jj*3+1:jj*3+3) = I(x-2+floor(jj/5),y-2+mod(jj,5),:);
                end                
            otherwise
        end
    end
    
    % �ϲ�rgbtmp������ɫѵ�����ݼ�
    rgbNotFish=[rgbNotFish;rgbtmp];
    
    prompt1 = {'go on?'};
    name1 = 'Segmenting manual';
    numlines1 = 1;
    defaultanswer1={'yes'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    in = inputdlg(prompt1,name1,numlines1,defaultanswer1,options);
    try
        go = in{1};
    catch
        go = 'no';
    end
end
save([outd,'rgbNotFish.mat'],'rgbNotFish');

%% Begin ANN Training
input = [rgbFish; rgbNotFish];
output = [ones(size(rgbFish,1),1); zeros(size(rgbNotFish,1),1)];
nnstart
