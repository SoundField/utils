clear;clc;
% --------------------------路径---------------------------------------
[fidd,filepath] = uigetfile('.line','open mat file');

%--------------------------TL曲线---------------------------------------
%%
tl=load([filepath,fidd]);                %修改到文件所在目录下 
figure,plot(tl(1:length(tl),1)/1000,tl(1:length(tl),2),'LineWidth',1.5)
axis ij;grid on;
r = tl(1:length(tl),1);
% hold on;plot(r/1000,10*log10(r),'k');
% hold on;plot(r/1000,15*log10(r),'k');
hold on;
legend('Rams','20log(r)');
xlabel('Range(km)','fontsize',16);
ylabel('Transmission Loss(dB)','fontsize',14);
set(gca,'FontSize',14);
saveas(gca,[filepath,fidd(1:end-5),'.fig']);