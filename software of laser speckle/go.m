clear all;
clc;
loadlibrary('spec2contr.dll')
a=imread('w0.bmp');
%disp ����1
K =rgb2gray(a);
K=spec2contr(double(K),5);
%disp ����2
V=contr2velo_improved(size(K,1)*size(K,2),K);
%disp ����3
imagesc(a);%��ʾԭʼɢ��ͼ
axis off;
title('ԭʼɢ��ͼ')
disp ɢ��ͼ
%set(gca,'position',[0.05 0.11 0.7750 0.8150]);
figure,imagesc(K);%��ʾ�ռ�ı�ͼ
axis off;
title('�ռ�ı�ͼ')
caxis([0.08 0.8]);
disp �ı�ͼ
%set(gca,'position',[0.05 0.11 0.7750 0.8150]);
figure,imagesc(V);%��ʾ����ٶ�ͼ
axis off;
title('����ٶ�ͼ')
caxis([50 200]);
disp �ٶ�ͼ
%set(gca,'position',[0.05 0.11 0.7750 0.8150]);
