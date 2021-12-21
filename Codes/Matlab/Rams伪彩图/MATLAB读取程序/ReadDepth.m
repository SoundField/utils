function [A] = ReadDepth(Zr)
%%%  读取分界面处的声压和位移。  %%%
%%%  建立地形文件 geo.txt。  %%%
file1 = 'geo.txt'; file2 = 'time_print.txt';file3 = 'weicai.txt';    %需要修改到文件所在目录下
%%
fid = fopen(file2);
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
dr = str2double(file_a(12:20));     %水平步长【dr*ndr】 
ndz = str2double(file_a(40:42));
zmplt = str2double(file_a(50:58));
fclose(fid);
%%
dixing = load(file1);dz = dzf;
aa = load(file3);
x = 0:dr:(range1)*dr;                 %距离
y = zeros(1,depth-1);
lengthx = length(x);lengthy = 4;
c = reshape(aa,lengthy,lengthx);   %伪彩图数据矩阵
nzr = ceil(Zr/dr);
A = c(nzr,:);
%---------------------------------------------------------------------
