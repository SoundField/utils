%%Ƶɢ���gif
clc;clf;close all;

%�������
load Hw;

%% ȷ���׷�ͼ����ʽ����ָ�����⣬������������ʽ
aa=size(pp);NN=aa(2);p=pp;fs=1000;
Num_depth = [30];              %����������
Num_range = 1:1:aa(2);                         %ֻ����һ�������
for ll=1:1:NN
    T_min=(Num_range(ll))*dr/1500-1;   %ʱ�䴰���ޣ���������ֱ�ﲨ�ĵ���ʱ��.*exp(-1i*2*pi*T_min(l)*Freq)
    aa=reshape( pp(Num_depth,Num_range(ll),:) ,1,length(Freq));
    HH(ll,:)=aa.*exp(-1i*2*pi*T_min*Freq);
end
H_si=[zeros(NN,Time*floor(freq_up)-1),HH,zeros(NN,Time*floor(N/2/Time-freq_down))];%numr�У�N�� HH��Ƶ�������е�������
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

%ȷ��ͼ���ڲɼ��Ĺ����а��������ἰ����
ax = gca;
ax.Units = 'pixels';
pos = ax.Position;
ti = ax.TightInset;
rect = [-ti(1),-ti(2),pos(3)+ti(1)+ti(3),pos(4)+ti(2)+ti(4)];
%��ָ���ķ�Χ�ڻ��ͼ���ļ�
frame = getframe(ax,rect);
im = frame2im(frame);

%����gif�ļ���ָ������ʽ��д����֡ͼ��
k = 1;
%��Ԫ�������ɼ�����ͼ�񣬷�����淴תͼ����
[I{k},map{k}]=rgb2ind(im,256);
imwrite(I{k},map{k,1},'mygif.gif','gif','Loopcount',Inf,'DelayTime',0.2);
k = k + 1;


%��ͼ���ɼ���gif��
for i = 2:100:NN
    sig = R_flip(i,:)';[sS,sF,sT] = spectrogram(sig,kaiser(256,4),200,128,fs);
    pcolor(sT,sF,abs(sS));shading interp;ylim([0 100]);
    xlabel('Time');ylabel('Freq');grid on;
    str = ['Disper of Pekris Waveguide [created by xl]'];
    title({str},'Interpreter','latex');
    %����gif
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
%��ͼ���෴��˳����д��gif��
% for i = (k-1):-1:1
%     imwrite(I{i},map{i},'mygif.gif','gif','WriteMode','append','DelayTime',0.1);
% end