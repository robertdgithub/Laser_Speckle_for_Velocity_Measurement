
function  speckle
clear all

figure('MenuBar','none','Name','speckle','NumberTitle','off','Position',[500, 200, 900, 600]);

uicontrol('Style','PushButton','String','非成像散斑','Position',[3,530,90,40],'CallBack',@BasicPressed);
uicontrol('Style','PushButton','String','成像散斑','Position',[3,430,90,40],'CallBack',@Pressed);
uicontrol('Style','PushButton','String','动态成像散斑','Position',[3,330,90,40],'CallBack',@dPressed);
uicontrol('Style','PushButton','String','动态成像散斑时间积分','Position',[3,290,120,40],'CallBack',@tPressed);
uicontrol('Style','PushButton','String','动态散斑图相关性','Position',[3,130,90,40],'CallBack',@convPressed);

k=1024;
A=rand(k);
B=rand(k);

A2=ones(k);
A2(600:1024,:)=0.6;
A2(300:599,:)=0.2;

function BasicPressed(~,~)
%非成像散斑
%表面（0,1）均匀分布，相位（0,1）均匀分布
A1=rand(k);
Aobj=A1.*exp(1i*2*pi*rand(1024));
fou=fft2(Aobj);
inten=abs(fou).^2;
% figure

subplot(2,1,1)
mesh(inten)
title('圆形复高斯散斑')
colormap gray
colorbar
a=inten(:);
subplot(2,1,2)
hist(a,50000)
xlabel('圆形复高斯散斑的光强分布')
end


 function Pressed(~,~)
%成像散斑
%表面高度阶梯状
A2=ones(k);
A2(600:1024,:)=0.6;
A2(300:599,:)=0.2;


Aobj=A2.*exp(1i*2*pi*rand(1024));
%random各种分布
fou=fft2(Aobj);

% lens pupil
se=strel('disk',256);
H = getnhood(se);
m=H.*H;
g=zeros(511,256);
m=[g,m,g];
g=zeros(511,1);
g=[m,g];
m=zeros(256,1024);
g=[m;g;m];
m=zeros(1,1024);
g=[g;m];


% ？？
fou=g.*fftshift(fou);

fou=fft2(fou);
inten=abs(fou).^2;

subplot(2,1,1)
mesh(inten)
title('反射强度分区的表面散斑图')
colormap gray
colorbar
%散斑强度取对数后的分布
a=inten(:);
subplot(2,1,2)
hist(a,50000)
xlabel('表面散斑图的光强分布')
 end


function dPressed(~,~)
%%动态散斑：k=2，n=10时的一张散斑图


z=sqrt(-2*log(B)).*cos(2*pi*A+pi/2/9);
f=cdf('unid',z,4);

Aobj=A2.*exp(1i*2*pi*f);

%random各种分布
fou=fft2(Aobj);

se=strel('disk',300);
H = getnhood(se);
m=H.*H;
g=zeros(599,212);
m=[g,m,g];
m=[m,zeros(599,1)];
g=zeros(212,1024);
g=[g;m;g];
g=[g;zeros(1,1024)];

fou=g.*fou;

fou=fft2(fou);
inten=abs(fou).^2;
subplot(1,2,1)
mesh(inten)
view(0,90)
title('10张动态散斑图中的第2张')
colormap gray
colorbar
% 散斑强度取对数后的分布
a=inten(:);
subplot(1,2,2)
hist(a,50000)
axis([0,.5e12,0,1000])
xlabel('散斑图的光强分布')
end

function tPressed(~,~)
%%动态散斑50张时间积分散斑图


se=strel('disk',300);
H = getnhood(se);
m=H.*H;
g=zeros(599,212);
m=[g,m,g];
m=[m,zeros(599,1)];
g=zeros(212,1024);
g=[g;m;g];
g=[g;zeros(1,1024)];

inten=zeros(1024);
for o=1:25
z=sqrt(-2*log(A)).*cos(2*pi*B+pi*(o-1)/2/49);
f=cdf('unif',z,0,4);
Aobj=A2.*exp(1i*2*pi*f);
fou=fft2(Aobj);
fou=g.*fou;
fou=fft2(fou);
inten1=abs(fou).^2;
inten=inten1+inten;
end

figure
mesh(inten/10)
view(0,90)
title('50张动态散斑图中前25张的时间积分')
colormap gray
colorbar
end

function convPressed(~,~)
%第二张动态散斑图与第七张动态散斑图之间的相关性
k=256;
A=rand(k);
B=rand(k);


A2=ones(k);
A2(150:256,:)=0.6;
A2(50:149,:)=0.2;
% A2=rand(n);


se=strel('disk',100);
H = getnhood(se);
m=H.*H;
g=zeros(199,28);
m=[g,m,g];
m=[m,zeros(199,1)];
g=zeros(28,256);
g=[g;m;g];
g=[g;zeros(1,256)];

% cor=zeros(1,100);
% 
% z=sqrt(-2*log(A)).*cos(2*pi*B);
% f=cdf('unif',z,0,4);
% Aobj=A2.*exp(1i*2*pi*f);
% fou=fft2(Aobj);
% fou=g.*fftshift(fou);
% fou=fft2(fou);
% inten=abs(fou).^2;
% 
% for l=1:50
% 
%     z1=sqrt(-2*log(A)).*cos(2*pi*B+pi*(l-1)/2/49);
%     f1=cdf('unif',z1,0,4);
%     Aobj1=A2.*exp(1i*2*pi*f1);
%     fou1=fft2(Aobj1);
%     fou1=g.*fftshift(fou1);
%     fou1=fft2(fou1);
%     inten1=abs(fou1).^2;
%     
%     cor(1,l)=corr2(inten,inten1);
% end
% plot(cor)  
% title('50张动态散斑图')
% ylabel('第1张与第K张散斑图的相关系数')
% xlabel('第K张散斑图')
% axis([0,50,0.6,1.05])

z=sqrt(-2*log(A)).*cos(2*pi*B+pi/2/49);
f=cdf('unif',z,0,4);
Aobj=A2.*exp(1i*2*pi*f);

%random各种分布


se=strel('disk',100);
H = getnhood(se);
m=H.*H;
g=zeros(199,28);
m=[g,m,g];
m=[m,zeros(199,1)];
g=zeros(28,256);
g=[g;m;g];
g=[g;zeros(1,256)];

fou=fft2(Aobj);
fou=g.*fftshift(fou);

fou=fft2(fou);
inten=abs(fou).^2;

z1=sqrt(-2*log(A)).*cos(2*pi*B+pi*2/2/49);
f1=cdf('unif',z1,0,4);
Aobj=A2.*exp(1i*2*pi*f1);

%random各种分布
fou=fft2(Aobj);
fou=g.*fftshift(fou);

fou=fft2(abs(fou).^2);
inten2=abs(fou).^2;
c1=xcorr2(inten,inten2);
figure
subplot(1,2,1)
mesh(c1*1e-20*3/2)
title('100张动态散斑图中第2、3张之间的相关性')
view(-90,0)
colorbar

c1=corr2(inten,inten2);


z3=sqrt(-2*log(A)).*cos(2*pi*B+pi*49/2/49);
f3=cdf('unif',z3,0,4);
Aobj=A2.*exp(1i*2*pi*f3);

%random各种分布
fou=fft2(Aobj);
fou=g.*fftshift(fou);

fou=fft2(abs(fou).^2);
inten3=abs(fou).^2;
c2=xcorr2(inten,inten3);
subplot(1,2,2)
mesh(c2*1e-20*7/5)
title('100张动态散斑图中第2、50张之间的相关性')
view(-90,0)
colorbar


end

end



