for k = 1
    k
BW = Bx{1,k};
binaryImage=BW;
% Skeletonize  
skeletonizedImage = bwmorph(binaryImage, 'skel', inf);  
% distance transform.  
Dist_Img = bwdist(~binaryImage);  
% multiply  
CenLine_Img = Dist_Img .* single(skeletonizedImage);  
% binarize  
T_level = 0.85;  
CenLine_Img = uint8(255*im2bw(CenLine_Img, T_level));  
  
% % display  
% figure(1),  
% imshow(skeletonizedImage, []);  
% figure(2),  
% imshow(Dist_Img, []);  
% figure(3),  
imshow(CenLine_Img);  
[a,b] = find(CenLine_Img~=0);
datato=[b,a];
%[idx,C,sumd,D] = kmedoids(datato,10);

[idx,C,sumd,D] = kmedoids(datato,10,'Distance','seuclidean');
%[idx,C,sumd,D] = kmedoids(xx,2);
imshow(BW);
hold on;
plot(C(:,1),C(:,2),'ro')
hold on;

plot(b,a,'b.')
% for i = 1:480
%   plot(B1{1,1}(i,1),B1{1,1}(i,2),'r.');
%   hold on;
%   pause(0.01);
%   plot(C(:,1),C(:,2),'bo')
% 
% end
pause(0.05);
hold off;
% hold on;
% plot(C(:,1),C(:,2),'bo')
end

 [L,num] = bwlabel(BW, 4);

   S0 = regionprops(L,'Area'); %求连通区域面积
    S0 = cat(1,S0.Area);
    C0 = regionprops(L,'Centroid'); %求连通区域位置
    C0 = cat(1,C0.Centroid);
    len0 = regionprops(L,'MajorAxisLength'); % 椭圆长轴
    len0 = cat(1,len0.MajorAxisLength);
    wid0 = regionprops(L,'MinorAxisLength'); % 椭圆短轴
    wid0 = cat(1,wid0.MinorAxisLength);
    Ec0 = regionprops(L,'Eccentricity'); % 离心率
    Ec0 = cat(1,Ec0.Eccentricity);
    B0=regionprops(L,'BoundingBox'); % 求连通区域边界
    B0=cat(1,B0.BoundingBox);
%     O0=regionprops(L,'Orientation'); % 椭圆方向
%     O0=cat(1,O0.Orientation);
    medianEc = median(Ec0); % 离心率中位数，平均值会被带偏
   