function [A] = ReadDepth(Zr)
%%%  ��ȡ�ֽ��洦����ѹ��λ�ơ�  %%%
%%%  ���������ļ� geo.txt��  %%%
file1 = 'geo.txt'; file2 = 'time_print.txt';file3 = 'weicai.txt';    %��Ҫ�޸ĵ��ļ�����Ŀ¼��
%%
fid = fopen(file2);
file_a = fgetl(fid);time_prog  = str2double( file_a(10:17) );    %��������ʱ��
file_a = fgetl(fid);time_start = str2double( file_a(9:17))*10^str2double( file_a(19:21) );  %����ʼʱ��
file_a = fgetl(fid);depth = str2double( file_a(15:18) );range1 = str2double( file_a(24:28) );   %��ȵ����;������
file_a = fgetl(fid);
dzw = str2double( file_a(8:16) );   %�Ǿ����������1   dzw
dzf = str2double( file_a(27:35) );  %�Ǿ����������2   dzf
dzb = str2double( file_a(51:59) );  %�Ǿ����������3   dzb
file_a = fgetl(fid);
H1 = str2double(file_a(8:15));      %����12�������
H2 = str2double(file_a(27:35));     %����23�������
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
fclose(fid);
%%
dixing = load(file1);dz = dzf;
aa = load(file3);
x = 0:dr:(range1)*dr;                 %����
y = zeros(1,depth-1);
lengthx = length(x);lengthy = 4;
c = reshape(aa,lengthy,lengthx);   %α��ͼ���ݾ���
nzr = ceil(Zr/dr);
A = c(nzr,:);
%---------------------------------------------------------------------
