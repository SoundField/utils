%  clear  all;
% close all;
% clc;

tl=load('E:\rams_0.5_Collins\tl.line');
tt=load('E:\rams_0.5_Collins\tl.line');
S=size(tl);
row=S(1,1);
column=S(1,2);
range=tl(4:row,1)/1000; 
tlc=tl(4:row,2);%-tt((4:200)-3,2);
figure;
plot(range,tlc,'b');
title(['z=' num2str(tl(3,1)) 'm transmission loss']);hold on
xlabel('r/km');
ylabel('TL/dB');
axis ij
grid on;
%  clear  all;
% close all;
% clc;
% 
% load('C:\Users\xcx\Desktop\rams_0.5_Collins\rams_0.5_Collins\tl.line');
% range=tl(4:20,1)/1000; 
% tlc=tl(4:20,2);
% figure;
% plot(range,tlc,'r');
% title('n=1');
% xlabel('r/km');
% ylabel('TL/dB');
% axis ij
% grid on;

% %%%mapping画图程序
% %%%%%%%%%TL曲线
% clear all
% %close all
% clc
% tl=load('C:\Users\xcx\Desktop\rams_0.5_Collins\rams_0.5_Collins\tl.line');                %修改到文件所在目录下 
% freq=tl(1,1);   Zs=tl(2,1);   Zr=tl(3,1);
% x=tl(4:length(tl),1)/1000;
% y=tl(4:length(tl),2);
% % x=downsample(x,10);
% % y=downsample(y,10);
% figure,plot(x,y)
% axis ij
% xlabel('Range(km)','fontsize',16);ylabel('Transmission Loss(dB)','fontsize',14);
% set(gca,'FontSize',14);
% legend('rams0.5')
% title(['f=' num2str(freq) 'Hz，Zs=' num2str(Zs) 'm，Zr=' num2str(Zr) 'm，p'],'fontsize',14);
% %  xlim([0,90])
% %%%%%%声场伪彩图
% load 'dt.txt';             %修改到文件所在目录下
% weicai=dt;
% freq=weicai(1,1);
% dr=weicai(2,1);              %一列数的第二个数为水平步长【dr*ndr】
% dr=dr/1000;
% ndr=weicai(3,1); 
% dz=weicai(4,1);                   %一列数的第三个数为垂直步长
% ndz=weicai(5,1);
% Zs=weicai(6,1);
% lz=weicai(9,1);                %深度上的点数
% aa=weicai(10:length(weicai),1);    %画伪彩的值
% range=length(aa)/lz;           %画伪彩时距离上的个数
% x=ndr*dr:ndr*dr:range*ndr*dr;                 %距离
% y=dz:ndz*dz:(lz-1)*ndz*dz+dz;              %深度
% lengthx=length(x);
% lengthy=length(y);
% c = reshape(aa,lengthy,lengthx);
% figure,pcolor(x,y,c)
% axis([dr*ndr range*dr*ndr dz 300]);        %设置伪彩图中水平和垂直距离
% axis ij
% shading interp
% jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb)
% caxis([20,120]);  %设置伪彩图中色棒范围
% colorbar('EastOutside')
% % title(['f=' num2str(freq) 'Hz，Zs=' num2str(Zs) 'm，p'],'fontsize',14);
% title(['f=' num2str(freq) 'Hz，Zs=' num2str(Zs) 'm  ， p'],'fontsize',20);
% xlabel('Range(m)','fontsize',22);
% ylabel('Depth(m)','fontsize',22);
% set(gca,'FontSize',14);
% 
% 
% % load 'w.txt';             %修改到文件所在目录下
% % weicai=w
% % freq=weicai(1,1);
% % dr=weicai(2,1);              %一列数的第二个数为水平步长【dr*ndr】
% % dr=dr/1000;
% % ndr=weicai(3,1); 
% % dz=weicai(4,1);                   %一列数的第三个数为垂直步长
% % ndz=weicai(5,1);
% % Zs=weicai(6,1);
% % lz=weicai(9,1);                %深度上的点数
% % aa=weicai(10:length(weicai),1);    %画伪彩的值
% % range=length(aa)/lz;           %画伪彩时距离上的个数
% % x=ndr*dr:ndr*dr:range*ndr*dr;                 %距离
% % y=dz:ndz*dz:(lz-1)*ndz*dz+dz;              %深度
% % lengthx=length(x);
% % lengthy=length(y);
% % c = reshape(aa,lengthy,lengthx);
% % figure,pcolor(x,y,c)
% % axis([dr*ndr range*dr*ndr dz 300]);        %设置伪彩图中水平和垂直距离
% % axis ij
% % shading interp
% % jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb)
% % caxis([20,120]);  %设置伪彩图中色棒范围
% % colorbar('EastOutside')
% % % title(['f=' num2str(freq) 'Hz，Zs=' num2str(Zs) 'm，p'],'fontsize',14);
% % title(['f=' num2str(freq) 'Hz，Zs=' num2str(Zs) 'm  ， p'],'fontsize',20);
% % xlabel('Range(km)','fontsize',22);
% % ylabel('Depth(m)','fontsize',22);
% % set(gca,'FontSize',14);
% % 
