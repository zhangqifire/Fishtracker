% ����ȥ��������Ⱥ���˶��켣׷���㷨
% ��ʹ�ÿ������˲�����ֱ������ʷƽ����Ԥ��
% ���¼����Ƶ��֡���ܹ����׼ȷ��
% �����Ƶ�ı����ȶ��������ڼ�СmaxError�Ӷ�����ڸ��ܶ��µ���ȷ��
% v1,2014.12.10,�ܹ�����ʶ��
% v2,2014.12.10,�޸������Ƶ֡��ʾΪʵ����Ƶ��Ӧ��֡��,�޸���ʶ�𵽵���ʵ����ٶ�ʧ�����
% v3,2014.12.10,�޸���ˮ��ָ������Cshape��fishsegn�е�bug
% v4,2014.12.11,�޸�v3��������ѭ��bug,����С�޸�,����maxJump
% v5,2014.12.11,�޸�����ָ�����۷���,�����txt�����Ӷ�Ӧ��Ƶ�е�֡λ��
% v6,2014.12.11,��������������Ƶ����,��������������resize,maxJump
% v7,2015.03.16,���Ӹ���Ȥ���򣬳�ʼλ���ֶ�ѡ��

% 2015.9.10 ��Ҫǡ������minLight��fishRatio���Բ�Ҫ��
% 2015.9.8  ��������Ž��ж��ηָ
            % �޸���Ծ������������֡��Խ�󣬿������������ΧԽ�󣻣�mycost.m��
            % �����Ѿ���ʧ�ĸ���ҲҪ�޶�������Χ���ڿ������������Χ�ڵ�costһ����֮���Inf ��mycost.m��
            % ���ڶ�ʧ�ĸ������³��ֵģ���������κ�λ���ǶԵģ�����������ɱ�����Ϊ0��ᵼ����ʶ�𵽵ĸ�������
            % ����ʵ�ʸ�������ʱ����ʧ�ĸ���ռ��δ��ʧ�����λ�ã���δ��ʧ�ĸ������µĶ�ʧ���塣
            % ��ˣ�A(ii,:)=9999999999(�㹻���������������Inf)�� ��mycost.m��
            % ������С��minLight�ĵ��Ϊȫ�ڣ�֮������graythresh������Ч���ܺã�fishsegn_s1.m��
% 2015.9.7  �޸ĻҶ�ת��ֵͼ��levelֵΪ0.05��fishsegn_s1.m��
            % ���ηָ��ͼ����Ȼ��С����С���һ�����ͨ����ȥ�� ��fishsegn_s2.m��
% 2015.9.6  ����Ϳ�˰�ɫ��Ҫʶ��ĸ������˸������֣����Ҹ����С��÷ǳ�С��
            % ������������϶�ʧ����������Ϊ��ֹ�ˡ���Ϊ�ܿ��������ϵ�λ�û���̬���⣬���¶�ʧ��
            % �����˶������о����ж�ʧҲ�Ƕ϶�������ֻ�о�ֹ�˲ſ���������֡��ʧ��
% 2015.6.11 �������Ƶ�����ö�������
% 2015.4.12 ��߶����Ŷ�ָ���ļ���׼ȷ��
% 2015.4.11 �޸�������ǰ�����Ȳ����Сʱ�Զ��Ҷ�ת��ֵͼ�����bug���޸���mycost�������Ȩ��
% 2015.4.10 �޸����Ŷ����һֱ��Conf0��bug
% 2015.04.09,�޸�ƥ������ΪLAP���⣬�������㷨+�Զ���ɱ�������������ԡ�����mainRgb.m�޸ġ�

clc; clear -regexp [^V,^region]; close all; fclose all;
% novideo = 0; % 1Ϊ���������Ƶ
testmode = 0; % 1Ϊ����ģʽ
testframe = Inf; % Ҫ��ʼ��ͣ�����֡
rgb = 0; % 0Ϊ����ȥ��������ʶ��
tic
disp('��ʼʶ��......')
disp('minLight����Ҫ����������ǡ������')
%##########��ʼ��������##########
% ��Ӧ�ø���ʵ���������mycost.m�и���Ȩ�ش�С��Ŀǰ2,1,2�����û��Ǻܺ����
v95=10; %��������������v99��95%����������ٶ�С��v95,��ʼ����д��㣬��40����������һ�Σ�����testVip.m�����ȡ������v95
minS=50; %��S1,��ȡ��ʽͬv95��һ��ʼ��������С�㣬��20��ʵ��ֻҪ���������ʶ��Ĳ��־͹���
maxS=400; %vmax��S99����ȡ��ʽͬv95��һ��ʼ�������ô�㣬��500������Ϊ�Դ�����ʵ�������������ʡ�
           %��ʶ�����Smax�������ص����ֵ������������ʵ�������������S95��S99��
% fishRatio = 0.01; % fishRatio����������ռ������������ı�����һ���Դ󣬵�Խ�ӽ�Խ�á�ͨ���ñ���С��0.01��
minLight = 10; % ������С��minLight�ĵ��Ϊȫ�ڣ�֮������graythresh������Ч���ܺã�fishsegn_s1.m��
                % ����ʵ����Ƶ�лҶ�ͼ�����������������ֵ��ȷ��minLight���ȱ������µĸ������ص��������ֵ��С����
con=0;
[ffile,fpath] = uigetfile({'*.avi';'*.mp4'},'Select the Video file');
prompt = {'fish number','start frame','frame number','frame interval',...
    'dark fish?','strelSize','resize','maxJump','isabort','novideo'};
name = 'Input Parameters';
numlines = 1;
% minS,maxSӦ��С�㣬��Ϊ��������������С���
defaultanswer={'20','1','0','1','yes','1','1','1','0','0'}; 
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
in = inputdlg(prompt,name,numlines,defaultanswer,options);
num0 = str2double(in{1}); % ��֪������
sf = str2double(in{2}); % ��ʼʶ���֡�����ڴ�֡����interval֡�в������ص�
N = str2double(in{3});  % ʶ���֡��,���Ϊ0��ʶ��������Ƶ
interval = str2double(in{4}); % ���ݴ���֡���
darkfish = in{5}; % ��ȱ�������Ϊyes������Ϊno
strelSize = str2double(in{6}); % ��̬ѧ�����ṹԪ�ش�С��һ�����ȡ4
resize = str2double(in{7}); % ���ű���
maxJump = str2double(in{8}); % �����Ծϵ�����������൱����Χ
isabort = str2double(in{9}); % Ϊ1��ʾ���ͬʱ������2����ΪNaN���������
novideo = str2double(in{10}); % 1������Ƶ������������Ƶ
%##########������������##########
fname=[fpath,ffile];
if testmode==1
    outd=[fpath,'output\test-',ffile,'\'];
else
    outd=[fpath,'output\',ffile,'\'];
end
set_minLight = 0; % 2015.9.10 ��һ��ʶ�����Ƶ����Ҫ����minLight
if(~isdir(outd)) 
    mkdir(outd);
    % 2015.9.10 ��һ��ʶ�����Ƶ����Ҫ����minLight
    set_minLight = 1;
end

% ��ȡͼ���ɰ�
% region=imread([fname,'_�ɰ�.png']);
% region=imresize(region,resize); %����



try
    if ~(exist('V','var') && strcmp(ffile,V.name))
        V=VideoReader(fname); % ��ȡ��Ƶǰ���ж���Ƶ�Ƿ��Ѿ�����ȡ�������ڵ���ʱ���ٵȴ���ȡ��Ƶ��ʱ��
    end
catch
    V=VideoReader(fname); 
end
try
    region=imread([fname,'_�ɰ�.png']);
    %region=rgb2gray(region);
catch exception   
 uiwait(msgbox('ѡ����Ȥ����,�����ɰ�','Success','modal'));
   
region = roipoly(read(V,100));
region =uint8(~ region.*255);
imwrite(region , [fname,'_�ɰ�.png']);
end
try
    %����Ѿ����ڱ���ͼƬ��ֱ�Ӷ�ȡ������ͼ��㱳��ͼƬ
    back=imread([fname,'_background_',num2str(resize),'.png']);
catch exception
    back=background(V,resize);
    imwrite(back,[fname,'_background_',num2str(resize),'.png']);
end

if N==0 || sf+(N-1)*interval > V.NumberOfFrames
    % ���NΪ0����Ҫʶ������֡��������Ƶ��֡������N����Ϊ���ֵ
    N = (V.NumberOfFrames - sf) / interval + 1;
end
% ��ʼ����������ļ�
fname=[outd,ffile]; %�ļ���ȥ��׺
fid=strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.txt');
fid2=strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_lost.txt'); %��֡��Ϣ
fid2=fopen(fid2,'w'); 

% д������ļ�
fid=fopen(fid,'w');
% % д�����' Confidence   px1    py1   px2   py2...'
fprintf(fid, ' FramePos  Confidence ');
for line=1:num0
    fprintf(fid,['   px',int2str(line),'      py',int2str(line),'   ']);
end
fprintf(fid,'\r\n');

if novideo ~= 1 % ����Ƶģʽ����ͼ��������Ƶ
% ��ʼ�������Ƶ
aviobj = VideoWriter(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N)),'MPEG-4');
aviobj.FrameRate = V.FrameRate;
if V.FrameRate>100
aviobj.FrameRate = 25;
end
open(aviobj);  
end

% ��ʼ��Ҫ��¼�ı���
px=NaN(num0,N); % x����λ��
py=px; % y����λ��
vx=px; % x�����ٶ�
vy=px; % y�����ٶ�
ax=px; % x������ٶ�
ay=py; % y������ٶ�
len=px; % ��Բ��ϵõ����峤
wid=px; % ��Բ��ϵõ������
Ec=px; % ��Բ������
S=px; % ������
Conf=NaN(1,N); % ׼ȷ������ָ��

countLost = zeros(num0,1); % ��ʼ����֡������
abort = 0; % ��ʼ����ֹ�����ʶ��
target_last = zeros(num0,1); % ��¼��һ����ƥ����
% ���������տ�ʼ�����֡
% ������ʱҪ�������ص�����Ҫ��֡���棬�������ֹ��ָ�
 frame=sf;
 save('test.mat');
 %����������С����ģ��������û��Ŀ��Ĵ���
 save('continue.mat','con');
 parameters
 while con==0
     pause(2);
     load continue.mat
 end
 load para.mat
 for frame=sf:interval:sf+(N-1)*interval
    colNow = ceil((frame-sf+1)/interval); %��ǰ֡�ڴ��������е����
    I0=read(V,frame); % ԭʼ��Ƶ֡ͼ��
    I0=rgb2gray(I0);
    I0=imresize(I0,resize); %����   
    % �Ҷ�ģʽ����
    if strcmp(darkfish,'yes')
        I = back - I0; %�뱳������Ϊ��������������back-I����֮I-back
    else
        I = I0 - back;
    end
    I=I-region; 
    
%%     % 2015.9.10 ������minLight���޸���fishsegn_s1����Ҫ�ⲽ��
%     if colNow == 1
%         % ����I�����ȸߵĲ������㣬���ȵ�(�ӽ�0��һ�㶼��С��10)�Ĳ����Ǳ���
%         % ���Ƕ��ڱ����������Ȳ���С����������ֱ���Զ�ת��Ϊ��ֵͼЧ�����ܻ�ܲ� 
%         % �����ڴ���ָ��ǰ������ǰ���ͱ����Ĳ��� % 2015.4.11
%         Ivalue = zeros(1,10);
%         Ivaluetmp = reshape(I,1,[]);
%         if fishRatio > 0.5 || fishRatio <= 0
%             disp('fishRatio���ô��󣬽���Ĭ��ֵ0.01���м���')
%             fishRatio = 0.01;
%         end
%         for tmpii = 1:length(Ivalue)
%             Ivalue(tmpii) = prctile(Ivaluetmp, 100*((1-fishRatio) + (tmpii-1)/length(Ivalue)*fishRatio));
%         end
%         Ivalue = mean(Ivalue);
%     end
%     % ���ڵ�Ivalue�����Ͽ�����Ϊǰ���ͱ������ȵķֽ��ߣ���Ivalue���Ȼ�С�Ŀ϶����Ǳ����ˡ�
%     I(I<Ivalue)=0; % �������߻Ҷ�ת��ֵͼ����ȷ��
    
    % 2015.9.10 ��һ��ʶ�����Ƶ����Ҫ����minLight
    if set_minLight == 1
        imshow(I)
        title('�ҳ�ǰ���ͱ����ķֽ�����ֵ�����潫������ֵ���ڸ���ֵ�ĵ�������Ϊ0������ߺ���ת��ֵͼ��׼ȷ��')
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
        manual; % 2015-3-17 �˹��ָ�ȷ����ʼ�������ص�
    end
    
    if testmode==1 && frame>=testframe
        keyboard
    end
   
    
%    imwrite(I,sprintf('fig/%d.png',frame));
    [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s1( I,num0,strelSize,minS,maxS,minLight ); %����ʶ����    
    if colNow > 3
        % �Ѿ����˳�ʼ���Ľ׶Σ����Ǹ���ƥ�������ηָ�
        % ��Ԥƥ��һ�Σ����ݽ���Է�Χ�ڴ��������ֵ�ļ��Ž��зָ�
        % ��ɱ����� source�У�target��     
        A = mycost(px(:,colNow-1),py(:,colNow-1),vx(:,colNow-2),vy(:,colNow-2),S(:,colNow-1), ...
            ax(:,colNow-3),ay(:,colNow-3),len(:,colNow-1),C0(:,1),C0(:,2),S0,v95,maxJump,countLost);
        % ��һ����ʶ��Ϊ������ûʲô���⣬����ʶ�𲻳����Ƚ��Ѱ죬����Ҫ����ʶ���������
        % ͬһ��������������֡�б���ʶ����ĸ��ʽ�С�����Կ�����һ������ʶ��
        target_indices = assignmentoptimal(A); % �˴�AΪ�ɱ���ԽСԽ��
        ind = find(target_indices==0); % ����ûƥ���ϵ�
        while ~isempty(ind)
            % �������û��ƥ���ϵģ���Դ��ڸø�����ܷ�Χ�ڵ���������зָ�            
            indtmp = ind(1); % ����indtmp����
            ind(1) = []; % ��ind(1)�Ƴ�
            % ��һ����λ��ΪԲ�ģ���v95+maxJump*lenΪ�뾶
            pxtmp = px(indtmp,colNow-1); %x���꣬ͼ���е�x,y��
            pytmp = py(indtmp,colNow-1); %y����
            distmp = v95 + maxJump * len(indtmp,colNow-1); % ���ܷ�Χ��С
            disC0 = sqrt((pxtmp-C0(:,1)).^2 + (pytmp-C0(:,2)).^2); % ��C0�������
            C0indtmp = find(disC0<=distmp); % ��Ҫ�ָ����ͨ�������
            if ~isempty(C0indtmp)
                % ��������Ҫ�ָ���Ųŷָ��ȻҪ����
                [Ec0,len0,wid0,S0,C0,Conf0,Iseg] = fishsegn_s2( Iseg,C0indtmp,Conf0,num0,minS );
            end
        end
    end
    Conf(colNow)=Conf0; % ׼ȷ������ָ��,�����������
    
    if colNow==1
        px(:,colNow)=C0(:,1); %x���꣬ͼ���е�x,y��
        py(:,colNow)=C0(:,2); %y����
        len(:,colNow)=len0; % ��Բ��ϵõ����峤
        wid(:,colNow)=wid0; % ��Բ��ϵõ������
        Ec(:,colNow)=Ec0; % ��Բ������
        S(:,colNow)=S0; % ������
        Conf(colNow)=Conf0; % ׼ȷ������ָ��
    else
        if colNow<=3
            source = [px(:,colNow-1),py(:,colNow-1)];
            target = [C0(:,1),C0(:,2)];
            % ����Ǹտ�ʼ������hungarianlinker.m����ƥ��
            % ��ʼ������������������������ƥ�䣬����Ҫ��ʶ��׼ȷ
            target_indices = hungarianlinker(source, target); 
        else
            % ��ɱ����� source�У�target��         
            A = mycost(px(:,colNow-1),py(:,colNow-1),vx(:,colNow-2),vy(:,colNow-2),S(:,colNow-1), ...
                ax(:,colNow-3),ay(:,colNow-3),len(:,colNow-1),C0(:,1),C0(:,2),S0,v95,maxJump,countLost);
              %%%%%%% test for A %%%%%%%%%%%%
%             frame
%             [A,Av,Aa,As] = mycost(px(:,colNow-1),py(:,colNow-1),vx(:,colNow-2),vy(:,colNow-2),S(:,colNow-1), ...
%                 ax(:,colNow-3),ay(:,colNow-3),len(:,colNow-1),C0(:,1),C0(:,2),S0,v95,maxJump,countLost);
            
            % ��һ����ʶ��Ϊ������ûʲô���⣬����ʶ�𲻳����Ƚ��Ѱ죬����Ҫ����ʶ���������
            % ͬһ��������������֡�б���ʶ����ĸ��ʽ�С�����Կ�����һ������ʶ��
            target_indices = assignmentoptimal(A); % �˴�AΪ�ɱ���ԽСԽ��
            
            %%%%% ���ݳɱ�����������Ŷ� %%%%%%%%%%%
            % ������ƥ��Ľ�����ϴβ�һ�����ҳɱ�������������������򽵵����Ŷ�
            if ~isequal(target_indices, target_last)
                Aminr = min(A,[],2); % A��ÿһ�е���Сֵ
                B = ( A - repmat(Aminr,1,size(A,2))) ./ repmat(Aminr,1,size(A,2)) ; % B�ǳɱ����������ƶ�
                if sum(sum(B<0.1,2) >= 2) >= 2
                    Conf(colNow)=Conf(colNow)-30; % �����2��������2�����ϵĿ���ƥ�䣬������׳���
                end
            end
            
            if testmode==1 && frame>=testframe
                disp(A)
                disp(target_indices')
            end
            target_last = target_indices; % ������һ����ƥ����
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % ����֡�������Ϣ����ʶ������        
        indyes = find(target_indices~=0); % ����ƥ���ϵ�
        while ~isempty(indyes)
            % ����û�ж�֡���㣬����������״̬
            indtmp = indyes(1); % ����indtmp����
            indyes(1) = []; % ��indyes(1)�Ƴ�
            countLost(indtmp) = 0; % ���ö�֡������
            px(indtmp,colNow)=C0(target_indices(indtmp),1); %x���꣬ͼ���е�x,y��
            py(indtmp,colNow)=C0(target_indices(indtmp),2); %y����
            len(indtmp,colNow)=len0(target_indices(indtmp)); % ��Բ��ϵõ����峤
            wid(indtmp,colNow)=wid0(target_indices(indtmp)); % ��Բ��ϵõ������
            Ec(indtmp,colNow)=Ec0(target_indices(indtmp)); % ��Բ������
            S(indtmp,colNow)=S0(target_indices(indtmp)); % ������
            vx(indtmp,colNow-1) = px(indtmp,colNow) - px(indtmp,colNow-1); % �ٶ�
            vy(indtmp,colNow-1) = py(indtmp,colNow) - py(indtmp,colNow-1);
            if colNow>2
                ax(indtmp,colNow-2) = px(indtmp,colNow) + px(indtmp,colNow-2) - 2 * px(indtmp,colNow-1); % ���ٶ�
                ay(indtmp,colNow-2) = py(indtmp,colNow) + py(indtmp,colNow-2) - 2 * py(indtmp,colNow-1); 
            end
        end
       
        ind = find(target_indices==0); % ����ûƥ���ϵ�
        while ~isempty(ind)
            % ��������㶪ʧ������ʧ�����Ϣ����ʷ��ϢԤ�����
            disp(['�����˶�֡��֡λ��Ϊ��',num2str(frame)])
            fprintf(fid2,['�����˶�֡��֡λ��Ϊ��',num2str(frame),'\r\n']);
            % ��������Ӧ���ٶȴ�С����
            vel = sqrt((vx(:,colNow - 2).^2 + vy(:,colNow -2).^2)); % t-2ʱ���ٶȴ�С
            vxtmp = vx(:,colNow-2)*2 - vx(:,colNow-3);
            vytmp = vy(:,colNow-2)*2 - vy(:,colNow-3);
            vtmp = sqrt(vxtmp.^2 + vytmp.^2); % t-1ʱ���ٶȴ�С
            pxtmp = px(:,colNow-1) + vxtmp .* vel ./ vtmp;
            pytmp = py(:,colNow-1) + vytmp .* vel ./ vtmp;
            
            indtmp = ind(1); % ����indtmp����
            ind(1) = []; % ��ind(1)�Ƴ�
            countLost(indtmp) = countLost(indtmp) + 1;
            if countLost(indtmp) > 2
                % ����Ԥ����2֡������ο�ʼдNaN,Ԥ��3֡ͨ��������
                Conf (colNow)=Conf(colNow)-20; % ׼ȷ������ָ��,ÿ��ʧһ�����Сһ��ֵ��������ʧ��֡Ҫ����������
                px(indtmp,colNow)=NaN; %x���꣬ͼ���е�x,y��
                py(indtmp,colNow)=NaN; %y����
                len(indtmp,colNow)=NaN; % ��Բ��ϵõ����峤
                wid(indtmp,colNow)=NaN; % ��Բ��ϵõ������
                Ec(indtmp,colNow)=NaN; % ��Բ������
                S(indtmp,colNow)=NaN; % ������
                vx(indtmp,colNow-1) = NaN; % �ٶ�
                vy(indtmp,colNow-1) = NaN;
                if colNow>2
                    ax(indtmp,colNow-2) = NaN; % ���ٶ�
                    ay(indtmp,colNow-2) = NaN; 
                end
            else
                Conf(colNow)=Conf(colNow)-10; % ׼ȷ������ָ��,ÿ��ʧһ�����Сһ��ֵ
                px(indtmp,colNow)=pxtmp(indtmp); %x���꣬ͼ���е�x,y��
                py(indtmp,colNow)=pytmp(indtmp); %y����
                len(indtmp,colNow)=len(indtmp,colNow-1); % ��Բ��ϵõ����峤
                wid(indtmp,colNow)=wid(indtmp,colNow-1); % ��Բ��ϵõ������
                Ec(indtmp,colNow)=Ec(indtmp,colNow-1); % ��Բ������
                S(indtmp,colNow)=S(indtmp,colNow-1); % ������
                vx(indtmp,colNow-1) = px(indtmp,colNow) - px(indtmp,colNow-1); % �ٶ�
                vy(indtmp,colNow-1) = py(indtmp,colNow) - py(indtmp,colNow-1);
                if colNow>2
                    ax(indtmp,colNow-2) = px(indtmp,colNow) + px(indtmp,colNow-2) - 2 * px(indtmp,colNow-1); % ���ٶ�
                    ay(indtmp,colNow-2) = py(indtmp,colNow) + py(indtmp,colNow-2) - 2 * py(indtmp,colNow-1); 
                end
            end
        end

        if sum(isnan(px(:,colNow)))>=2 && isabort==1
            % �����2�����������Ѿ����׶�ʧ����ֹ����
            abort = 1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    end
    
    % % д��ʶ����
    fprintf(fid,'   %5.1f       %5.1f    ',frame,Conf(colNow));
    for line=1:num0
        fprintf(fid,'%7.2f  %7.2f  ',px(line,colNow),py(line,colNow));
    end
    fprintf(fid,'\r\n');
    
    if novideo ~= 1 % ����Ƶģʽ����ͼ��������Ƶ
    if testmode==1 % ����ģʽ
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
        linehistory = 20; %���ߵĳ���
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
    if abort==1
        break
    end
end
%% �������ս����MAT�ļ�
px(:,colNow+1:end)=[]; % ����δʶ�𲿷֣���С�ļ���С
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

tElapsed = toc % ��ʾ�����������������ʱ��
fprintf(fid2,['tElapsed = ',num2str(tElapsed),' s\r\n']);
% ���������MAT�ļ�
save(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_par.mat'), ...
    'num0','sf','N','interval','darkfish','V',...
    'strelSize','tElapsed','resize','maxJump','isabort','v95','minS','maxS','minLight');
if novideo ~= 1 % ����Ƶģʽ����ͼ��������Ƶ
close(aviobj);
end
close('all');
fclose('all');
se = frame; %ʵ��ʶ������һ֡
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.txt'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'.txt'))
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.mat'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'.mat'))
if novideo ~= 1 % ����Ƶģʽ����ͼ��������Ƶ
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'.mp4'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'.mp4'))
end
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_lost.txt'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'_lost.txt'))
movefile(strcat(fname,'_sf=',num2str(sf),'_N=',num2str(N),'_par.mat'),...
    strcat(fname,'_sf=',num2str(sf),'_ef=',num2str(frame),'_par.mat'))
load handel.mat; %���������������ʾ��
sound(y,2*Fs);

disp('��������')