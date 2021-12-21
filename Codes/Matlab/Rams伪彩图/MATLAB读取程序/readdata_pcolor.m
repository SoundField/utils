%----------------------  均匀网格  -------------------------------------
clear;close all;clc;
% --------------------------路径---------------------------------------
[fidd,filepath] = uigetfile('.txt','open txt file');
% fidd='weicai.txt';
% filepath='C:\Users\51158\Desktop\论文\程序\Rams1.5\Rams1.5\Rams1.5\';
file1 = 'geo.txt';file2 = 'time_print.txt';file3 = fidd;
% -----------------------读取声场信息-------------------------------------
fid = fopen([filepath,file2]);
file_a = fgetl(fid);time_prog  = str2double( file_a(9:20) );    %程序运行时间
file_a = fgetl(fid);time_start = str2double( file_a(9:17))*10^str2double( file_a(19:21) );  %程序开始时间
file_a = fgetl(fid);depth = str2double( file_a(15:18) );range1 = str2double( file_a(24:28) );   %深度点数和距离点数
file_a = fgetl(fid);
dz = str2double( file_a(8:16) );   %均匀网格分区1   dzw
file_a = fgetl(fid);
freq = str2double(file_a(15:22));
zs   = str2double(file_a(44:54));
file_a = fgetl(fid);
zr   = str2double(file_a(18:25));
rmax = str2double( file_a(46:53) ); %最大计算距离
file_a = fgetl(fid); 
dr = str2double(file_a(12:20));     %水平步长【dr*ndr】 
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
%% -----------------------声场伪彩图------------------------------------
weicai = load([filepath,fidd]);             %修改到文件所在目录下
depth=weicai(1,1);                %一列数的第一个数为深度上的个数
dr=weicai(2,1)/1000;                   %一列数的第二个数为水平步长【dr*ndr】
dz=weicai(3,1);                   %一列数的第三个数为垂直步长
zmplt=weicai(4,1);                %一列数的第四个数为画伪彩图时取的深度
aa=weicai(5:length(weicai),1);    %画伪彩的值
range = length(aa)/(depth+1);
x=dr:dr:range*dr;                 %距离
y=0:dz:(depth)*dz;
lengthx=length(x);
lengthy=length(y);
c = reshape(aa,lengthy,lengthx);
figure,pcolor(x,y,c)
axis([dr range*dr 0 zmplt-100]);        %设置伪彩图中水平和垂直距离
axis ij;shading interp
jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb);
caxis([30,120]);  %设置伪彩图中色棒范围
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
ii=1;       %伪彩矩阵列数
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