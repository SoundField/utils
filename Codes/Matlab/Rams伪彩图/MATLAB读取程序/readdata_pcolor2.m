%----------------------  非均匀网格  ------------------------------------
clear;close all;clc;
% -----------------------------------------------------------------------
[fidd,filepath] = uigetfile('.txt','open txt file');
%%
file1 = 'geo.txt';file2 = 'time_print.txt';file3 = fidd;
% -----------------------读取声场信息-------------------------------------
fid = fopen([filepath,file2]);
file_a = fgetl(fid);time_prog  = str2double( file_a(10:17) );    %程序运行时间
file_a = fgetl(fid);time_start = str2double( file_a(9:17))*10^str2double( file_a(19:21) );  %程序开始时间
file_a = fgetl(fid);depth = str2double( file_a(15:18) );range1 = str2double( file_a(24:28) );   %深度点数和距离点数
file_a = fgetl(fid);
dzw = str2double( file_a(8:16) );   %非均匀网格分区1   dzw
dzf = str2double( file_a(27:35) );  %非均匀网格分区2   dzf
dzb = str2double( file_a(51:59) );  %非均匀网格分区3   dzb
file_a = fgetl(fid);
H1 = str2double(file_a(8:15));      %分区12划分深度
H2 = str2double(file_a(27:35));     %分区23划分深度
file_a = fgetl(fid);
freq = str2double(file_a(15:22));
zs   = str2double(file_a(44:54));
file_a = fgetl(fid);
zr   = str2double(file_a(18:25));
rmax = str2double( file_a(46:53) ); %最大计算距离
file_a = fgetl(fid); 
dr = str2double(file_a(12:20))/1000;     %水平步长【dr*ndr】 
ndz = str2double(file_a(40:42));
zmplt = str2double(file_a(50:58));
fclose(fid);
%% -----------------------声场伪彩图------------------------------------
weicai = load([filepath,fidd]);             %修改到文件所在目录下
aa=weicai(1:length(weicai),1);    %画伪彩的值
range1 = rmax/1000/dr;
x=dr:dr:range1*dr;                 %距离
for i=1:1:ceil(H1/dzw)
    h(i)=dzw;
end
for i=ceil(H1/dzw)+1:1:ceil((H2-H1)/dzf)+ceil(H1/dzw)
    h(i)=dzf;
end
for i=ceil((H2-H1)/dzf)+ceil(H1/dzw)+1:1:depth-1
    h(i)=dzb;
end
for i=2:1:depth-1
    y1(1)=0;
    y1(i)=y1(i-1)+h(i-1);
end
y = y1(1:depth-1);
lengthx=length(x);
lengthy=length(y);
c = reshape(aa,lengthy,lengthx);
figure,pcolor(x,y,c)
axis([dr range1*dr 0 zmplt]);        %设置伪彩图中水平和垂直距离
axis ij;shading interp
jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb);
caxis([30,120]);  %设置伪彩图中色棒范围
colorbar('NorthOutside')
xlabel('Range(km)','fontsize',16);
ylabel('Depth(m)','fontsize',14);
set(gca,'FontSize',14);
saveas(gca,[filepath,'weicai.png']);
% close all;