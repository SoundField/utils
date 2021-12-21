close all;clear;clc;warning off;tic;
filepath = 'C:\Users\51158\Desktop\����\����\����RAMs�������������ʱ��';
%------------------------------------------------------------------------
%%%% ����Ӵ��������������ź�Ƶ��
f=50;w=2*pi*f;              % �����ź�Ƶ��
fs=1000;                  % ����Ƶ��
N=1000;     n=0:(N-1);      Time=N/fs;% ��t=1/fs,��f=1/Time    
T1=0:1/fs:4/f;
T2=4/f+1/fs:1/fs:Time;
T=[T1 T2];
Sig=[1/2*sin(w*T1).*(1-cos(w*T1/4)),zeros(1,length(T2))];
figure;         %% ��ʱ����
plot(T,Sig);xlim([0 Time]);

%%%%% ����CW���ź�Ƶ��
% f=30;w=2*pi*f;              % �����ź�Ƶ��
% fs=1000;                    % ����Ƶ��
% N=15000;     n=0:(N-1);      Time=N/fs;% ��t=1/fs,��f=1/Time  
% T=0:1/fs:Time-1/fs;Sig=zeros(1,length(T));      % ʱ�䳤��5s
% cnt=100;                      % ������ڸ���
% t=0:1/fs:cnt/f-1/fs;
% sig=sin(w*t); % ����hanning����ȥ�԰�
% Sig(1:length(sig))=sig;
% figure;         %% ��ʱ����
% plot(T,Sig);
% xlabel('Time��s��');ylabel('Pressure��Pa��');
% ylim([-1.1 1.1]);
%% fft
xlabel('ʱ�䣨s��');ylabel('���ȣ�V��');
xlim([0 T(end)]);ylim([-1.0 1.0]);
Sig_p=fft(Sig,N);      %% ��Ƶ��
l=-(N/2):1:N/2-1;F=l/Time;
plot(F,fftshift(abs((Sig_p))));
%% ------------example-1:Munk��������---------------------------
freq = 50; w = 2*pi*freq; zs = 50;zr = 50;
dr = 10.0; ndr = 1; rmax = 5000;
c0 = 1500.0; np = 6; irot = 0; theta = 60.0;
zmax = 1500; dz = 1.0;ndz = 1; zmplt = 1500;

%%%%%%%%%low acoustic velocity elastic layer(������)λ����Ϣ
lavel_depth=[0.0    100.0
5000.0  100.0];
%%%%%%%%% SSP
ssp = [0.00     1500.0];
%%%%  cp1 cs1 rhob1 attnp1 attns1  %%%%
cp = [0.0  1600.0];        
cs = [0.0 800.0];
rhob = [0.0 1.5]; 
attnp = [0.0 0.1
500.0   1.0
1000.0  10.0]; 
attns = [0.0 0.2
500.0   2.0
1000.0  20.0];
csch=cst(freq,lavel_depth(1,2),cp(end),cs(end),rhob(end)*1000);
% disp(Cscholte);
%%%%%%%%%%%%%%%%%%%%% .in�ļ���д %%%%%%%%%%%%%%%%%%%%%%%%%%%%        
fid = fopen( 'rams.in', 'w' );
fprintf( fid, '%s \r\n', 'rams_S3.in' );
fprintf( fid, '%0.1f    %0.1f    %0.1f\r\n ',freq ,zs ,zr  );
fprintf( fid, '%0.1f  %0.1f  %i\r\n',rmax ,dr ,ndr);
fprintf( fid, '%0.1f  %0.1f  %i  %0.1f \r\n',zmax, dz, ndz, zmplt);
fprintf( fid, '%0.1f  %i  %i  %0.1f \r\n',c0, np, irot, theta );
%-------------------------------------------------------------------------
fprintf( fid, '%0.1f       %0.1f\r\n', lavel_depth(1,:)' );
fprintf( fid, '%0.1f       %0.1f\r\n', lavel_depth(2:end,:)');
fprintf( fid, '%i %i\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f\r\n', ssp(1,:)' );
fprintf( fid, '%0.1f       %0.1f\r\n', ssp(2:end,:)');
fprintf( fid, '%i %i\r\n',-1,-1);
% -------------------------------------------------------------------
fprintf( fid, '%0.1f       %0.1f\r\n', cp(1,:)' );
fprintf( fid, '%0.1f       %0.1f\r\n', cp(2:end,:)');
fprintf( fid, '%i %i\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f\r\n', cs(1,:)' );
fprintf( fid, '%0.1f       %0.01f\r\n', cs(2:end,:)');
fprintf( fid, '%i %i\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f\r\n', rhob(1,:)' );
fprintf( fid, '%0.1f       %0.01f\r\n', rhob(2:end,:)');
fprintf( fid, '%i %i\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f\r\n', attnp(1,:)' );
fprintf( fid, '%0.1f       %0.01f\r\n', attnp(2:end,:)');
fprintf( fid, '%i %i\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f\r\n', attns(1,:)' );
fprintf( fid, '%0.1f       %0.01f\r\n', attns(2:end,:)');
fprintf( fid, '%i %i\r\n',-1,-1);
%----------------------------------------------------------------------
%--------------------- .in�ļ���д��� --------------------------------%
fclose( fid );
%Time_axes=T_min+n/fs;       %ʱ������һ��������ʱ�䴰
%Freq_axes=l*2*pi*fs/N;      %Ƶ������˫�ߵĽ�Ƶ��   �����㺣����ѧӢ�İ�P.614��
freq_up=25;freq_down=75;    %����ŵ�������
Freq=freq_up:1/Time:freq_down;   %Ƶ�ʼ��1
for ii=1:length(Freq)
    waitbar(ii/length(Freq));
    freq=Freq(ii);
    fid = fopen( 'rams.in','r+' );
    fseek(fid,10,'bof');
    fprintf( fid, '\r\n %0.1f   ', freq );    % �ı�rams�����ļ���freq    
    fclose( fid );

    % %%%%%%%% ����rams����       
 %   !rams-����.exe
    system([filepath,'\','Rams_tj.exe']);
    aa = load('pcolor.txt');
    bb = load('weicai_pr.txt');
    cc = load('weicai_pi.txt');
    [p] = Readfjm(lavel_depth,aa,bb,cc);
    pp(:,ii)=p';     %size(p)Ϊ �������*Ƶ����  1000*328(40*8.19)
    dd = load('tl_ri.line');
    prs(:,ii) = dd(:,2)+1i*dd(:,3); 
end;
%%
numr = rmax/dr;Hw = prs;numf = (freq_down-freq_up)*Time+1;
H_si=[zeros(numr,freq_up*Time-1),Hw,zeros(numr,N/2-freq_down*Time)];%numr�У�1000��
H_se = conj(seqreverse(H_si));H=[H_si H_se];
RF=zeros(size(H));
clear aa bb cc dd p;
%% IFFT tmin = range/c0
ff=zeros(1,numf);
k=1;
for p=freq_up:1/Time:freq_down
    ff(k)=p;
    k=k+1;
end
fff=[zeros(1,24),ff,zeros(1,925)];
tt=(0:N-1)/fs;
for i=1:floor(numr/5):numr   %������پ�����٣������˷�ʱ��
    tmin=i*dr/1500+0.5;
    RF(i,:)=Sig_p.*H(i,:).*exp(-1i*tmin*2*pi*fff);
    R(i,:)=ifft(RF(i,:),N);
    ma=max(max(R(i,:)));
    R_f(i,:)=R(i,:)/ma;
end 
figure
for i=1:floor(numr/5):numr
    plot(tt,real(R_f(i,:))+3*i/floor(numr/10),'k');
    hold on
end
title('���20m��ͬ����Ĳ���')
xlabel('�Ա�ʱ��t-r/1493(s)')
ylabel('���루km��')


%% ifft(���ٸ���Ҷ���任)
% numr=rmax/dr;
% % p=conj(p);
% Num=[500,1500,2000,2500,3000];NN=length(Num);HH=zeros(NN,length(Freq));
% for l=1:1:NN
%     T_min(l)=(Num(l))*dr/1500-1;   %ʱ�䴰���ޣ���������ֱ�ﲨ�ĵ���ʱ��.*exp(-1i*2*pi*T_min(l)*Freq)
%     HH(l,:)=pp(Num(l),:).*exp(-1i*2*pi*T_min(l)*Freq);
% end
% 
% H_si=[zeros(NN,Time*floor(freq_up)-1),HH,zeros(NN,Time*floor(N/2/Time-freq_down))];%numr�У�N�� HH��Ƶ�������е�������
% H_se=conj(seqreverse(H_si));H=[H_si H_se];
% CC=cst(freq,H1,cp1(:,2),cs1(:,2),rhob1(:,2)*1000);  %Scholte������
% omega_l=(0:1:N-1)*2*pi/Time;
%��j��ʱ��� j��t��Tmin��Tmin+(N-1)/fs;
%��f=1/Time,f=[(-N/2)/Time,1/Time,(N/2-1)/Time];omega_l=f*2*pi;
%ÿ��ʱ����0��(N/2-1)�������ã�ÿ�����p�������о���㴦����Ƶ�㡣
% R=zeros(size(H));R_flip=zeros(size(H));RF=zeros(size(H));
% for i=1:NN   %i=1 2 3 4 5 6
% %     waitbar(i/NN);
%     RF(i,:)=Sig_p.*H(i,:);
%     R(i,:)=ifft(RF(i,:),N);
%     ma=max(max(abs(R(i,:))));
%     R(i,:)=R(i,:)/ma;
%     for m=1:N
%     R_flip(i,m)=R(i,N-m+1);
%     end
% end 

%%%
% nn=2;
% RRR=ifft(RF(nn,:),N);
%     for m=1:N
%     Rr_flip(m)=RRR(N-m+1);
%     end
% plot(n/fs,real(Rr_flip));
%%%


% figure
% for i=1:NN
%     waitbar(i/NN);
%     plot(n/fs,real(1*R_flip(i,:))+1.5*i);
%     hold on
% end
% grid on;
% title(['���',num2str(zr),'m��ͬ����Ĳ���'])
% text(39,32,[' F  =  ',num2str(f),' Hz',sprintf('\n'),'SD = ',num2str(zs),' m',sprintf('\n'),'RD = ',num2str(zr),' m'],'horiz','center','fontsize',7.5)
% xlabel('�Ա�ʱ��t-r/1500+1(s)')
% ylabel('���루km��');xlim([0 10]);
%% STFT(��ʱ����Ҷ�任)
% %%%ʱ��T {0��1/fs��Time}    Ƶ��f_stft {freq_up-10:1/Time:freq_down+10}
% %%%��ѹ p(���룬Ƶ��) example��7000x101
% Num=[500,1000,1500,2000,2500,3000];NN=length(Num);HH=zeros(NN,length(Freq));
% for l=1:1:NN
%     T_min(l)=(Num(l))*dr/1500-1;   %ʱ�䴰���ޣ���������ֱ�ﲨ�ĵ���ʱ��.*exp(-1i*2*pi*T_min(l)*Freq)
%     HH(l,:)=p(Num(l),:).*exp(-1i*2*pi*T_min(l)*Freq);
%     
% end
% H_si=[zeros(NN,Time*floor(freq_up)-1),HH,zeros(NN,Time*floor(N/2/Time-freq_down))];%numr�У�N�� HH��Ƶ�������е�������
% H_se=conj(seqreverse(H_si));H=[H_si H_se];
% R=zeros(size(H));R_flip=zeros(size(H));RF=zeros(size(H));
% for i=1:NN   %i=1 2 3 4 5 6
% %     waitbar(i/NN);
%     RF(i,:)=Sig_p.*H(i,:)*fs*(N-1)/N;
%     R(i,:)=ifft(RF(i,:),N);
%     ma=max(max(abs(R(i,:))));
%     R(i,:)=R(i,:)/ma;
%     for m=1:N
%     R_flip(i,m)=R(i,N-m+1);
%     end
% end 
% figure (3)
% nn=5;
% sig = R_flip(nn,:)';[sS,sF,sT] = spectrogram(sig,kaiser(256,4),200,128,fs);
% pcolor(sT,sF,abs(sS));shading interp;ylim([0 100]);xlim([0.2 9]);
% xlabel('Reduce Time/ t-r/1500+1');ylabel('Frequency/Hz');
% title([num2str(Num(nn)*dr/1000),'km���봦���ʱƵͼ']);


