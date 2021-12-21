%%频散输出gif
clc;clf;close all;

%获得数据
load Hw;

%% 确定首幅图的样式，并指定标题，坐标轴标题等样式
aa=size(pp);NN=aa(2);p=pp;fs=1000;
Num_depth = [30];              %多个接收深度
Num_range = 1:1:aa(2);                         %只能有一个距离点
for ll=1:1:NN
    T_min=(Num_range(ll))*dr/1500-1;   %时间窗下限，理论上是直达波的到达时间.*exp(-1i*2*pi*T_min(l)*Freq)
    aa=reshape( pp(Num_depth,Num_range(ll),:) ,1,length(Freq));
    HH(ll,:)=aa.*exp(-1i*2*pi*T_min*Freq);
end
H_si=[zeros(NN,Time*floor(freq_up)-1),HH,zeros(NN,Time*floor(N/2/Time-freq_down))];%numr行，N列 HH从频率上限列到下限列
H_se=conj(seqreverse(H_si));H=[H_si H_se];
R=zeros(size(H));R_flip=zeros(size(H));RF=zeros(size(H));
for i=1:NN
%     waitbar(i/NN);
    RF(i,:)=Sig_p.*H(i,:)*fs*(N-1)/N;
    R(i,:)=ifft(RF(i,:),N);
    ma=max(max(abs(R(i,:))));
    R(i,:)=R(i,:)/ma;
    for m=1:N
    R_flip(i,m)=R(i,N-m+1);
    end
end 
sig = R_flip(1,:)';[sS,sF,sT] = spectrogram(sig,kaiser(256,4),200,128,fs);
pcolor(sT,sF,abs(sS));shading interp;ylim([0 100]);
xlabel('Time');ylabel('Freq');grid on;

[a,b] = size(sS);
x = sT; y = sF' ;t = T;x1 = zeros(a,b);y1 = zeros(a,b);
str = ['Disper of Pekris Waveguide [created by xl]'];
title({str},'Interpreter','latex');

%确保图像在采集的过程中包括坐标轴及标题
ax = gca;
ax.Units = 'pixels';
pos = ax.Position;
ti = ax.TightInset;
rect = [-ti(1),-ti(2),pos(3)+ti(1)+ti(3),pos(4)+ti(2)+ti(4)];
%在指定的范围内获得图像文件
frame = getframe(ax,rect);
im = frame2im(frame);

%创建gif文件，指定其样式，写入首帧图像
k = 1;
%用元胞存贮采集到的图像，方便后面反转图像用
[I{k},map{k}]=rgb2ind(im,256);
imwrite(I{k},map{k,1},'mygif.gif','gif','Loopcount',Inf,'DelayTime',0.2);
k = k + 1;


%画图并采集到gif中
for i = 2:100:NN
    sig = R_flip(i,:)';[sS,sF,sT] = spectrogram(sig,kaiser(256,4),200,128,fs);
    pcolor(sT,sF,abs(sS));shading interp;ylim([0 100]);
    xlabel('Time');ylabel('Freq');grid on;
    str = ['Disper of Pekris Waveguide [created by xl]'];
    title({str},'Interpreter','latex');
    %制作gif
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1),-ti(2),pos(3)+ti(1)+ti(3),pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect);
    im = frame2im(frame);
    [I{k},map{k}]=rgb2ind(im,256);
    %
    imwrite(I{k},map{k},'mygif.gif','gif','WriteMode','append','DelayTime',0.1);
    k=k+1;
end
%将图像按相反的顺序再写到gif中
% for i = (k-1):-1:1
%     imwrite(I{i},map{i},'mygif.gif','gif','WriteMode','append','DelayTime',0.1);
% end