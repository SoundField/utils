

close all;
clear;
clc;
warning off;
tic;
%%%%% 发射声信号频率
f=50;w=2*pi*f;              % 发射信号频率
fs=1000;                  % 采样频率
N=10000;     n=0:(N-1);      Time=N/fs;% Δt=1/fs,Δf=1/Time    
T1=0:1/fs:4/f;T2=4/f+1/fs:1/fs:(N-1)/fs;T=[T1 T2];
Sig=[1/2*sin(w*T1).*(1-cos(w*T1/4)),zeros(1,length(T2))];
figure;         %% 画时域波形
plot(T,Sig);
%% fft
xlabel('时间（s）');ylabel('幅度（V）');
xlim([0 T(end)]);ylim([-1.0 1.0]);
Sig_p=fft(Sig,N);      %% 画频谱
l=-(N/2):1:N/2-1;F=l/Time;
figure;
plot(F,fftshift(abs((Sig_p))));
%%
%%%%%%%%%example-1:Munk环境参数
freq = 25; w = 2*pi*freq; zs = 85;zr = 99.0;
dr = 5.0; ndr = 1; rmax = 35000;H1 = 90.0; H2 = 110.0;
c0 = 1500.0; np = 8; irot = 0; theta = 60.0;
zmax = 1000; dzw = 1.0; dzf = 1.0; dzb = 1.0; ndz = 1; zmplt = 1000;

%%%%%%%%%low acoustic velocity elastic layer(沉积层)位置信息
lavel_depth=[0.00 100.00
            35000 100.00];

%%%%%%%%%high acoustic velocity elastic layer(弹性层)位置信息
havel_depth=[0.00 105.00
            35000 105.00];
%%%%%%%%% SSP
ssp = [0.00     1500.0
        100.0      1500.00];
%%%%  cp1 cs1 rhob1 attnp1 attns1  %%%%
cp1 = [0.0  1600.0];        
cs1 = [0.0 800.0];
rhob1 = [0.0 1.5]; 
attnp1 = [0.0 0.1
    100.0 0.1]; 
attns1 = [0.0 0.2
    100.0 0.2];
%%%%  cp2 cs2 rhob2 attnp2 attns2  %%%%
cp2 = [0.0  1600.0];        
cs2 = [0.0 800.0];
rhob2 = [0.0 1.5]; 
attnp2 = [100.0 0.1
    1000.0 10.0]; 
attns2 = [100.0 0.2
    1000.0 20.0];

% Cscholte=cst(freq,H1,cp1(end),cs1(end)，rhob1(end));
% disp(Cscholte);
%%%%%%%%%%%%%%%%%%%%% .in文件编写 %%%%%%%%%%%%%%%%%%%%%%%%%%%%        
fid = fopen( 'rams_S3.in', 'w' );
fprintf( fid, '%s \r\n', 'rams_S3.in' );
fprintf( fid, '%0.1f    %0.1f    %0.1f ',freq ,zs ,zr  );
fprintf( fid, '                          %s  %s  %s\r\n','freq' ,'zs' ,'zr');
fprintf( fid, '%0.1f  %0.1f  %i  %0.1f  %0.1f  %s  %s  %s  %s  %s\r\n'...
    ,rmax ,dr ,ndr ,H1 ,H2 ,'rmax' ,'dr' ,'ndr' ,'H1' ,'H2');
fprintf( fid, '%0.1f  %0.1f  %0.1f  %0.1f  %i  %0.1f',...
    zmax, dzw, dzf, dzb, ndz, zmplt);
fprintf( fid, '   %s  %s  %s  %s  %s  %s\r\n', 'zmax', 'dzw', 'dzf',... 
   'dzb', 'ndz', 'zmplt');
fprintf( fid, '%0.1f  %i  %i  %0.1f', c0, np, irot, theta );
fprintf( fid, '     %s  %s  %s  %s\r\n\r\n', 'c0',...
    'np', 'irot', 'theta' );
%%%%%%%%%%%%%%%%%%%%%
fprintf( fid,'%s \r\n',...
    '*******************************************************************');
fprintf( fid,'%s \r\n','*(1) layer:2      沉积层      bathymetry_1');
fprintf( fid,'%s \r\n\r\n',...
    '*******************************************************************');
fprintf( fid, '%0.1f       %0.1f', lavel_depth(1,:)' );
fprintf( fid, '      %s\r\n', 'bathymetry_1');
fprintf( fid, '%0.1f       %0.1f\r\n\r\n', lavel_depth(2:end,:)');
fprintf( fid, '%i     %i\r\n',-1,-1);
%%%%%%%%%%%%%%%%%%%%%%
fprintf( fid,'%s \r\n',...
    '*******************************************************************');
fprintf( fid,'%s \r\n','*(2) layer:3      弹性层     bathymetry_2');
fprintf( fid,'%s \r\n\r\n',...
    '*******************************************************************');
fprintf( fid, '%0.1f       %0.1f', havel_depth(1,:)' );
fprintf( fid, '      %s\r\n', 'bathymetry_2');
fprintf( fid, '%0.1f       %0.1f\r\n\r\n', havel_depth(2:end,:)');
fprintf( fid, '%i     %i\r\n',-1,-1);
%%%%%%%%%%%%%%%%%%%%%%
fprintf( fid,'%s \r\n',...
    '*******************************************************************');
fprintf( fid,'%s \r\n','*(3) layer:1      Water        cw  ');
fprintf( fid,'%s \r\n\r\n',...
    '*******************************************************************');
fprintf( fid, '%0.1f       %0.1f', ssp(1,:)' );
fprintf( fid, '      %s\r\n', 'z cw');
fprintf( fid, '%0.1f       %0.1f\r\n\r\n', ssp(2:end,:)');
fprintf( fid, '%i      %i\r\n',-1,-1);
%%%%%%%%%%%%%%%%%%%%%%%
fprintf( fid,'%s \r\n',...
    '*******************************************************************');
fprintf( fid,'%s \r\n','*(4) layer:2     cp1, cs1, rhob1, attnp1, attns1');
fprintf( fid,'%s \r\n\r\n',...
    '*******************************************************************');
fprintf( fid, '%0.1f       %0.1f', cp1(1,:)' );
fprintf( fid, '      %s\r\n', 'cp1');
fprintf( fid, '%0.1f       %0.1f\r\n', cp1(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', cs1(1,:)' );
fprintf( fid, '      %s\r\n', 'cs1');
fprintf( fid, '%0.1f       %0.01f\r\n', cs1(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', rhob1(1,:)' );
fprintf( fid, '      %s\r\n', 'rhob1');
fprintf( fid, '%0.1f       %0.01f\r\n', rhob1(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', attnp1(1,:)' );
fprintf( fid, '      %s\r\n', 'attnp1');
fprintf( fid, '%0.1f       %0.01f\r\n', attnp1(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', attns1(1,:)' );
fprintf( fid, '      %s\r\n', 'attns1');
fprintf( fid, '%0.1f       %0.01f\r\n', attns1(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
%%%%%%%%%%%%%%%%%%%%%%%
fprintf( fid,'%s \r\n',...
    '*******************************************************************');
fprintf( fid,'%s \r\n','*(5) layer:3     cp2, cs2, rhob2, attnp2, attns2');
fprintf( fid,'%s \r\n\r\n',...
    '*******************************************************************');
fprintf( fid, '%0.1f       %0.1f', cp2(1,:)' );
fprintf( fid, '      %s\r\n', 'cp1');
fprintf( fid, '%0.1f       %0.1f\r\n', cp2(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', cs2(1,:)' );
fprintf( fid, '      %s\r\n', 'cs1');
fprintf( fid, '%0.1f       %0.01f\r\n', cs2(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', rhob2(1,:)' );
fprintf( fid, '      %s\r\n', 'rhob1');
fprintf( fid, '%0.1f       %0.01f\r\n', rhob2(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', attnp2(1,:)' );
fprintf( fid, '      %s\r\n', 'attnp1');
fprintf( fid, '%0.1f       %0.01f\r\n', attnp2(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
fprintf( fid, '%0.1f       %0.1f', attns2(1,:)' );
fprintf( fid, '      %s\r\n', 'attns1');
fprintf( fid, '%0.1f       %0.01f\r\n', attns2(2:end,:)');
fprintf( fid, '%i     %i\r\n\r\n',-1,-1);
%%%%%%%%%%%%%%%%%%%%% .in文件编写完成 %%%%%%%%%%%%%%%%%%%%%%%%%%
fclose( fid );
T_min=sqrt((rmax)^2+(zr-zs)^2)/max(max(cp1(1,2),cs1(1,2)),c0);
freq_up=25;freq_down=70;
Freq=freq_up:1:freq_down;
% Sig_p_r=Sig_p;
p=zeros(rmax/dr,length(Freq));
for ii=1:length(Freq)
    waitbar(ii/length(Freq));
    freq=Freq(ii);
    fid = fopen( 'rams_S3.in','r+' );
    fseek(fid,10,'bof');
    fprintf( fid, '\r\n %0.1f   ', freq );    % 改变rams输入文件中freq    
    fclose( fid );

    % %%%%%%%% 调用rams计算       
    !rams.exe
    
    load('p_zr.txt');
    p(:,ii)=p_zr(:,2)+1i*p_zr(:,3);     %size(p)为 距离点数*频点数  1000*328(40*8.19)
    
end;

%% ifft
numr=rmax/dr;
% p=conj(p);
omega_l=linspace(freq_up,freq_down,length(Freq));
Num=[500,1000,1500,2000,2500,3000];NN=length(Num);
for l=1:1:NN
    T_min = (Num(l))*dr/1492-0.2;
    HH(l,:)=p(Num(l),:).*exp(-1i*2*pi*omega_l*T_min);
end
H=[zeros(NN,floor(freq_up-1)),HH,zeros(NN,floor(N-freq_down))];%numr行，1000列 从16列到46列是不同频率的

CC=cst(freq,H1,cp1(:,2),cs1(:,2),rhob1(:,2)*1000);
omega_l=(0:1:N-1)*2*pi/Time;
%第j个时间点 jΔt从Tmin到Tmin+(N-1)/fs;
%Δf=1/Time,f=[(-N/2)/Time,1/Time,(N/2-1)/Time];omega_l=f*2*pi;
%每个时间由0到(N/2-1)的求和求得，每个求和p都是所有距离点处单个频点。
R=zeros(size(H));R_flip=zeros(size(H));RF=zeros(size(H));
for i=1:NN   %i=1 101 201...901
    waitbar(i/NN);
    RF(i,:)=Sig_p.*H(i,:);
    R(i,:)=ifft(RF(i,:),N);
    ma=max(max(abs(R(i,:))));
    R(i,:)=R(i,:)/ma;
    for m=1:N
    R_flip(i,m)=R(i,N-m+1);
    end
end 

figure
for i=1:NN
    waitbar(i/NN);
    plot(n/fs,2*real(R_flip(i,:))+5*i);
    hold on
end
grid on;
title(['深度',num2str(zr),'m不同距离的波形'])
text(16,43,[' F  =  ',num2str(f),' Hz',sprintf('\n'),'SD = ',num2str(zs),' m',sprintf('\n'),'RD = ',num2str(zr),' m'],'horiz','center','fontsize',7.5)
xlabel('对比时间t-r/1500(s)')
ylabel('距离（km）')
% saveas(gca,['D:\非均匀网格\宽带Scholte波','\zs=',num2str(zs),'深度',num2str(zr),'m不同距离的波形','.fig']);
% figure
% plot(n/fs,(R_flip(2,:)));hold on;
T_pro(zs)=toc;
% saveas(gca,['D:\非均匀网格\宽带Scholte波','\zs=',num2str(zs),' ',num2str(rmax),'处接收的时域波形','.fig']);

%%
% % %%%%%%%%%%  ifft到时域 %%%%%%%
%  Sig_p_ri = seqreverse(Sig_p_r);
%  Sig_p_ri = conj(Sig_p_ri);              %%%%%%%%%%%%%%%取共扼
%  Sig_p_r = [Sig_p_r Sig_p_ri];    
%  sig_p_r=ifft(Sig_p_r);
% %%%%%%%% 截取信号补前 %%%%%%%%%
%  sig_p_r=[sig_p_r(6000:10000) sig_p_r(1:5999)];
%  baoluo1=hilbert(sig_p_r/max(abs(sig_p_r)));
%  figure;
%  plot(T,sig_p_r/max(abs(sig_p_r)));hold on;
% %  plot(T,abs(baoluo),'r')
%  title('声压');
% xlim([6.5 8]);
% 
% % Sig_vu_ri = seqreverse(Sig_vu_r);
% % Sig_vu_ri = conj(Sig_vu_ri);              %%%%%%%%%%%%%%%取共扼
% % Sig_vu_r = [Sig_vu_r Sig_vu_ri];    
% % sig_vu_r=ifft(Sig_vu_r);
% % % sig_vu_r=[sig_vu_r(6000:10000) sig_vu_r(1:5999)];
% % baoluo2=hilbert(sig_vu_r/max(abs(sig_vu_r)));
% % % figure;
% % % plot(T,sig_vu_r/max(abs(sig_vu_r)));hold on;
% % % plot(T,abs(baoluo),'r')
% % % title('水平振速');
% % % xlim([6.5 8]);
% % 
% % Sig_vv_ri = seqreverse(Sig_vv_r);
% % Sig_vv_ri = conj(Sig_vv_ri);              %%%%%%%%%%%%%%%取共扼
% % Sig_vv_r = [Sig_vv_r Sig_vv_ri];    
% % sig_vv_r=ifft(Sig_vv_r);
% % % sig_vv_r=[sig_vv_r(6000:10000) sig_vv_r(1:5999)];
% % baoluo3=hilbert(sig_vv_r/max(abs(sig_vv_r)));
% % % figure;
% % % plot(T,sig_vv_r/max(abs(sig_vv_r)));hold on;
% % % plot(T,abs(baoluo),'r')
% % % title('垂直振速');
% % % xlim([6.5 8]);
% % 
% % 
% T_dir=sqrt((zr-zs)^2+rmax^2)/c0;
% CC=cst(freq,H1,cp1(:,2),cs1(:,2));
% T_sch=(H1-zs)/c0+rmax/CC+(H1-zr)/c0;
% % 
% figure;
% subplot(1,1,1);
% plot(T,sig_p_r/max(abs(sig_p_r)));hold on;
% plot(T,abs(baoluo1),'r');
% xlim([0 20]);
% ylim([-1 1]);
% ylabel('归一化幅度');
% xlabel('时域到达时间t/s')
% % 
% T_pro=toc;
% % saveas(gca,['D:\非均匀网格\宽带Scholte波\','Zs=',num2str(zs),' 计算时间t=',num2str(T_pro),'.fig']);
% 
% % subplot(3,1,2);
% % plot(T+0.6,sig_vu_r/max(abs(sig_vu_r)));hold on;
% % plot(T+0.6,abs(baoluo2),'r');
% % xlim([7 7.6]);
% % ylim([-0.8 0.8]);
% % ylabel('归一化幅度');
% % 
% % subplot(3,1,3);
% % plot(T+0.6,sig_vv_r/max(abs(sig_vv_r)));hold on;
% % plot(T+0.6,abs(baoluo3),'r');
% % xlim([7 7.6]);
% % ylim([-0.8 0.8]);
% % ylabel('归一化幅度');