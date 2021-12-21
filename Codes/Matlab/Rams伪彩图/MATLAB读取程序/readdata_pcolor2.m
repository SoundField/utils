%----------------------  �Ǿ�������  ------------------------------------
clear;close all;clc;
% -----------------------------------------------------------------------
[fidd,filepath] = uigetfile('.txt','open txt file');
%%
file1 = 'geo.txt';file2 = 'time_print.txt';file3 = fidd;
% -----------------------��ȡ������Ϣ-------------------------------------
fid = fopen([filepath,file2]);
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
dr = str2double(file_a(12:20))/1000;     %ˮƽ������dr*ndr�� 
ndz = str2double(file_a(40:42));
zmplt = str2double(file_a(50:58));
fclose(fid);
%% -----------------------����α��ͼ------------------------------------
weicai = load([filepath,fidd]);             %�޸ĵ��ļ�����Ŀ¼��
aa=weicai(1:length(weicai),1);    %��α�ʵ�ֵ
range1 = rmax/1000/dr;
x=dr:dr:range1*dr;                 %����
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
axis([dr range1*dr 0 zmplt]);        %����α��ͼ��ˮƽ�ʹ�ֱ����
axis ij;shading interp
jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb);
caxis([30,120]);  %����α��ͼ��ɫ����Χ
colorbar('NorthOutside')
xlabel('Range(km)','fontsize',16);
ylabel('Depth(m)','fontsize',14);
set(gca,'FontSize',14);
saveas(gca,[filepath,'weicai.png']);
% close all;