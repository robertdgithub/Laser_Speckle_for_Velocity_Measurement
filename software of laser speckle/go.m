clear all;
clc;
loadlibrary('spec2contr.dll')
a=imread('w0.bmp');
%disp 调试1
K =rgb2gray(a);
K=spec2contr(double(K),5);
%disp 调试2
V=contr2velo_improved(size(K,1)*size(K,2),K);
%disp 调试3
imagesc(a);%显示原始散斑图
axis off;
title('原始散斑图')
disp 散斑图
%set(gca,'position',[0.05 0.11 0.7750 0.8150]);
figure,imagesc(K);%显示空间衬比图
axis off;
title('空间衬比图')
caxis([0.08 0.8]);
disp 衬比图
%set(gca,'position',[0.05 0.11 0.7750 0.8150]);
figure,imagesc(V);%显示相对速度图
axis off;
title('相对速度图')
caxis([50 200]);
disp 速度图
%set(gca,'position',[0.05 0.11 0.7750 0.8150]);
