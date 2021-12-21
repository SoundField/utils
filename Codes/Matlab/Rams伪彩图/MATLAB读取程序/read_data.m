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

% %%%mapping��ͼ����
% %%%%%%%%%TL����
% clear all
% %close all
% clc
% tl=load('C:\Users\xcx\Desktop\rams_0.5_Collins\rams_0.5_Collins\tl.line');                %�޸ĵ��ļ�����Ŀ¼�� 
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
% title(['f=' num2str(freq) 'Hz��Zs=' num2str(Zs) 'm��Zr=' num2str(Zr) 'm��p'],'fontsize',14);
% %  xlim([0,90])
% %%%%%%����α��ͼ
% load 'dt.txt';             %�޸ĵ��ļ�����Ŀ¼��
% weicai=dt;
% freq=weicai(1,1);
% dr=weicai(2,1);              %һ�����ĵڶ�����Ϊˮƽ������dr*ndr��
% dr=dr/1000;
% ndr=weicai(3,1); 
% dz=weicai(4,1);                   %һ�����ĵ�������Ϊ��ֱ����
% ndz=weicai(5,1);
% Zs=weicai(6,1);
% lz=weicai(9,1);                %����ϵĵ���
% aa=weicai(10:length(weicai),1);    %��α�ʵ�ֵ
% range=length(aa)/lz;           %��α��ʱ�����ϵĸ���
% x=ndr*dr:ndr*dr:range*ndr*dr;                 %����
% y=dz:ndz*dz:(lz-1)*ndz*dz+dz;              %���
% lengthx=length(x);
% lengthy=length(y);
% c = reshape(aa,lengthy,lengthx);
% figure,pcolor(x,y,c)
% axis([dr*ndr range*dr*ndr dz 300]);        %����α��ͼ��ˮƽ�ʹ�ֱ����
% axis ij
% shading interp
% jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb)
% caxis([20,120]);  %����α��ͼ��ɫ����Χ
% colorbar('EastOutside')
% % title(['f=' num2str(freq) 'Hz��Zs=' num2str(Zs) 'm��p'],'fontsize',14);
% title(['f=' num2str(freq) 'Hz��Zs=' num2str(Zs) 'm  �� p'],'fontsize',20);
% xlabel('Range(m)','fontsize',22);
% ylabel('Depth(m)','fontsize',22);
% set(gca,'FontSize',14);
% 
% 
% % load 'w.txt';             %�޸ĵ��ļ�����Ŀ¼��
% % weicai=w
% % freq=weicai(1,1);
% % dr=weicai(2,1);              %һ�����ĵڶ�����Ϊˮƽ������dr*ndr��
% % dr=dr/1000;
% % ndr=weicai(3,1); 
% % dz=weicai(4,1);                   %һ�����ĵ�������Ϊ��ֱ����
% % ndz=weicai(5,1);
% % Zs=weicai(6,1);
% % lz=weicai(9,1);                %����ϵĵ���
% % aa=weicai(10:length(weicai),1);    %��α�ʵ�ֵ
% % range=length(aa)/lz;           %��α��ʱ�����ϵĸ���
% % x=ndr*dr:ndr*dr:range*ndr*dr;                 %����
% % y=dz:ndz*dz:(lz-1)*ndz*dz+dz;              %���
% % lengthx=length(x);
% % lengthy=length(y);
% % c = reshape(aa,lengthy,lengthx);
% % figure,pcolor(x,y,c)
% % axis([dr*ndr range*dr*ndr dz 300]);        %����α��ͼ��ˮƽ�ʹ�ֱ����
% % axis ij
% % shading interp
% % jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb)
% % caxis([20,120]);  %����α��ͼ��ɫ����Χ
% % colorbar('EastOutside')
% % % title(['f=' num2str(freq) 'Hz��Zs=' num2str(Zs) 'm��p'],'fontsize',14);
% % title(['f=' num2str(freq) 'Hz��Zs=' num2str(Zs) 'm  �� p'],'fontsize',20);
% % xlabel('Range(km)','fontsize',22);
% % ylabel('Depth(m)','fontsize',22);
% % set(gca,'FontSize',14);
% % 
