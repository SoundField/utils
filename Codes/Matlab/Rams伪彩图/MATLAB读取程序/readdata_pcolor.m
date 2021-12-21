%----------------------  ��������  -------------------------------------
clear;close all;clc;
% --------------------------·��---------------------------------------
[fidd,filepath] = uigetfile('.txt','open txt file');
% fidd='weicai.txt';
% filepath='C:\Users\51158\Desktop\����\����\Rams1.5\Rams1.5\Rams1.5\';
file1 = 'geo.txt';file2 = 'time_print.txt';file3 = fidd;
% -----------------------��ȡ������Ϣ-------------------------------------
fid = fopen([filepath,file2]);
file_a = fgetl(fid);time_prog  = str2double( file_a(9:20) );    %��������ʱ��
file_a = fgetl(fid);time_start = str2double( file_a(9:17))*10^str2double( file_a(19:21) );  %����ʼʱ��
file_a = fgetl(fid);depth = str2double( file_a(15:18) );range1 = str2double( file_a(24:28) );   %��ȵ����;������
file_a = fgetl(fid);
dz = str2double( file_a(8:16) );   %�����������1   dzw
file_a = fgetl(fid);
freq = str2double(file_a(15:22));
zs   = str2double(file_a(44:54));
file_a = fgetl(fid);
zr   = str2double(file_a(18:25));
rmax = str2double( file_a(46:53) ); %���������
file_a = fgetl(fid); 
dr = str2double(file_a(12:20));     %ˮƽ������dr*ndr�� 
ndz = str2double(file_a(40:42));
zmplt = str2double(file_a(50:58));
file_a = fgetl(fid); 
c0 = str2double(file_a(7:15));
file_a = fgetl(fid); 
if file_a ~= -1
theta1 = str2double(file_a(10:20));
theta2 = str2double(file_a(34:42));
end
fclose(fid);
%% -----------------------����α��ͼ------------------------------------
weicai = load([filepath,fidd]);             %�޸ĵ��ļ�����Ŀ¼��
depth=weicai(1,1);                %һ�����ĵ�һ����Ϊ����ϵĸ���
dr=weicai(2,1)/1000;                   %һ�����ĵڶ�����Ϊˮƽ������dr*ndr��
dz=weicai(3,1);                   %һ�����ĵ�������Ϊ��ֱ����
zmplt=weicai(4,1);                %һ�����ĵ��ĸ���Ϊ��α��ͼʱȡ�����
aa=weicai(5:length(weicai),1);    %��α�ʵ�ֵ
range = length(aa)/(depth+1);
x=dr:dr:range*dr;                 %����
y=0:dz:(depth)*dz;
lengthx=length(x);
lengthy=length(y);
c = reshape(aa,lengthy,lengthx);
figure,pcolor(x,y,c)
axis([dr range*dr 0 zmplt-100]);        %����α��ͼ��ˮƽ�ʹ�ֱ����
axis ij;shading interp
jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb);
caxis([30,120]);  %����α��ͼ��ɫ����Χ
colorbar('NorthOutside')
xlabel('Range(km)','fontsize',16);
ylabel('Depth(m)','fontsize',14);
set(gca,'FontSize',14);
% saveas(gca,[filepath,'weicai.png']);
% close all;
%% gausess bea
if file_a ~= -1
for i=1:1:zmplt/dz
   z=dz:dz:zmplt;
   k0= 2*pi*freq/1500;
   p(i) = sqrt(k0)*tan(theta1)*exp(-k0^2/2*(z(i)-zs)^2*tan(theta1)^2)...
       *exp(1i*k0*(z(i)-zs)*sin(theta2));
   W2(i) = real(p(i))^2+imag(p(i))^2;
end
figure(2);
suptitle('Acoustic Pressure of Gaussian beaming ')
subplot(1,2,1);
plot(real(p),z,'LineWidth',0.1);
axis ij;ylabel('Depth,m');xlabel('Range = 0');
mp=max(p);md=z(find(p == mp));
text(mp/2,md+100*dz,[num2str(md),'m'])

subplot(1,2,2);
ii=1;       %α�ʾ�������
p_1=10.^( -c(:,ii)./20)*sqrt(dr*ii);
p_1( find(p_1==1e+20) )=min(p_1);
p_1=p_1/max(p_1);
plot(p_1,z,'LineWidth',0.1)
axis ij;xlabel('Range = dr');

figure(3);
suptitle('Angle of Gaussian beaming ')
% subplot(1,2,1);

% subplot(1,2,2);

sita=linspace(-45,45,length(p));
plot(sita,W2);
% mip=max(ip);mang=sita(find(ip == mip));
% text(mang*3/4,mip,[num2str(mang),'\circ']);
xlim([-45 45]);xlabel(['Angle,','\circ']);
end