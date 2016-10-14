%%%Prewitt Operator Corner Detection.m
%%%ʱ���Ż�--����������ȡ��ķ�����Harris�ǵ�
 %% 

%  Image = imread('884.jpg');                 % ��ȡͼ��
%  Image = im2uint8(rgb2gray(Image));   
  
Image=rgb2gray(I1);
dx = [-1 0 1;-1 0 1;-1 0 1];  %dx������Prewitt���ģ��
Ix2 = filter2(dx,Image).^2;   
Iy2 = filter2(dx',Image).^2;                                        
Ixy = filter2(dx,Image).*filter2(dx',Image);
 
%���� 9*9��˹���ڡ�����Խ��̽�⵽�Ľǵ�Խ�١�
h= fspecial('gaussian',11,2);
A = filter2(h,Ix2);       % �ø�˹���ڲ��Ix2�õ�A 
B = filter2(h,Iy2);                                 
C = filter2(h,Ixy);                                  
nrow = size(Image,1);                            
ncol = size(Image,2);                             
Corner = zeros(nrow,ncol); %zeros��������һ��ȫ����󣬹ʾ���Corner���������ѡ�ǵ�λ��,��ֵȫ�㣬ֵΪ1�ĵ��ǽǵ�

 %����t:��(i,j)������ġ����ƶȡ�������ֻ�����ĵ������������˸��������ֵ֮����
 %��-t,+t��֮�䣬��ȷ������Ϊ���Ƶ㣬���Ƶ㲻�ں�ѡ�ǵ�֮��
 t=20;
 
 %�Ҳ�û��ȫ�����ͼ��ÿ���㣬���ǳ�ȥ�˱߽���boundary�����أ�Ҳ���Ǵӵ�8�е�8�п�ʼ������
 %��Ϊ���Ǹ���Ȥ�Ľǵ㲢�������ڱ߽���
 %���˾�����һ�����ǵ���ҪĿ�����ҳ������ǽǵ�ĵ㣬��С��Χ���ӿ������ٶȡ�
 %����˼����������ĵ㣨i,j����Χ8��������7��8����Ҷ�ֵ��֮���ƣ���ô�����ĵ�Ӧ�ô���ƽ̹���򣬲�����Ϊ�ǵ㣬
 %������ĵ㣨i,j����Χֻ��1�������û�е���֮���ƣ���ô�����ĵ�Ҳ������Ϊ�ǵ㡣
 boundary=8;
 for i=boundary:nrow-boundary+1 
    for j=boundary:ncol-boundary+1
         nlike=0; %���Ƶ����
         if Image(i-1,j-1)>Image(i,j)-t && Image(i-1,j-1)<Image(i,j)+t 
            nlike=nlike+1;
         end
         if Image(i-1,j)>Image(i,j)-t && Image(i-1,j)<Image(i,j)+t  
            nlike=nlike+1;
         end
         if Image(i-1,j+1)>Image(i,j)-t && Image(i-1,j+1)<Image(i,j)+t  
            nlike=nlike+1;
         end  
        if Image(i,j-1)>Image(i,j)-t && Image(i,j-1)<Image(i,j)+t  
            nlike=nlike+1;
         end
         if Image(i,j+1)>Image(i,j)-t && Image(i,j+1)<Image(i,j)+t  
            nlike=nlike+1;
         end
         if Image(i+1,j-1)>Image(i,j)-t && Image(i+1,j-1)<Image(i,j)+t  
            nlike=nlike+1;
         end
         if Image(i+1,j)>Image(i,j)-t && Image(i+1,j)<Image(i,j)+t  
            nlike=nlike+1;
         end
         if Image(i+1,j+1)>Image(i,j)-t && Image(i+1,j+1)<Image(i,j)+t  
            nlike=nlike+1;
         end
         if nlike>=2 && nlike<=6
             Corner(i,j)=1;%�����Χ��2~6�����Ƶ㣬��(i,j)���ǽǵ�
         end;
     end;
 end;
CRF = zeros(nrow,ncol);    % CRF��������ǵ���Ӧ����ֵ,��ֵȫ��
 CRFmax = 0;                % ͼ���нǵ���Ӧ���������ֵ������ֵ֮�� 
k=0.05;   
% ����CRF
 %�����ϳ���CRF(i,j) =det(M)/trace(M)����CRF����ô��ʱӦ�ý������105�е�
 %����ϵ��k���ô�һЩ��k=0.1�Բɼ����⼸��ͼ����˵��һ���ȽϺ���ľ���ֵ
 for i = boundary:nrow-boundary+1
     for j = boundary:ncol-boundary+1
     if Corner(i,j)==1  %ֻ��ע��ѡ��
         M = [A(i,j) C(i,j);
              C(i,j) B(i,j)];      
         CRF(i,j) = det(M)-k*(trace(M))^2;
         if CRF(i,j) > CRFmax 
            CRFmax = CRF(i,j);    
        end;            
    end
 end;             
end;  
%CRFmax
 count = 0;       % ������¼�ǵ�ĸ���
 t=0.01;        
% ����ͨ��һ��3*3�Ĵ������жϵ�ǰλ���Ƿ�Ϊ�ǵ�
 for i = boundary:nrow-boundary+1 
for j = boundary:ncol-boundary+1
         if Corner(i,j)==1  %ֻ��ע��ѡ��İ�����
             if CRF(i,j) > t*CRFmax && CRF(i,j) >CRF(i-1,j-1) ......%?????ΪʲôҪCRF(i,j) > t*CRFmax����������֪
                && CRF(i,j) > CRF(i-1,j) && CRF(i,j) > CRF(i-1,j+1) ......
                && CRF(i,j) > CRF(i,j-1) && CRF(i,j) > CRF(i,j+1) ......
                && CRF(i,j) > CRF(i+1,j-1) && CRF(i,j) > CRF(i+1,j)......
                && CRF(i,j) > CRF(i+1,j+1) 
            count=count+1;%����ǽǵ㣬count��1
             else % �����ǰλ�ã�i,j�����ǽǵ㣬����Corner(i,j)��ɾ���Ըú�ѡ�ǵ�ļ�¼
                 Corner(i,j) = 0;     
            end;
         end; 
end; 
end; 
 %disp('�ǵ����');
 % disp(count)
 figure,imshow(Image);      % display Intensity Image
 hold on; 
% toc(t1)
 for i=boundary:nrow-boundary+1 
for j=boundary:ncol-boundary+1
         column_ave=0;
         row_ave=0;
         k=0;
        if Corner(i,j)==1
             for x=i-3:i+3  %7*7����
                 for y=j-3:j+3
                     if Corner(x,y)==1
 % ������ƽ������Ϊ�ǵ����꣬������ü���ƽ��������ƽ�����꣬�Խǵ����ȡ���岻��
                         row_ave=row_ave+x;
                         column_ave=column_ave+y;
                         k=k+1;
                     end
                 end
             end
         end
         if k>0 %��Χ��ֹһ���ǵ�           
             plot( column_ave/k,row_ave/k ,'g.');
         end
 end; 
end; 

