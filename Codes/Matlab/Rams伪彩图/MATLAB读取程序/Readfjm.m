function [p] = Readfjm(Depth,aa,bb,cc)
%%%  ��ȡ�ֽ��洦����ѹ��λ�ơ�  %%%
%%%  ���������ļ� geo.txt��  %%%
file2 = 'time_print.txt';    %��Ҫ�޸ĵ��ļ�����Ŀ¼��
%%
fid = fopen(file2);
file_a = fgetl(fid);time_prog  = str2double( file_a(10:17) );    %��������ʱ��
file_a = fgetl(fid);time_start = str2double( file_a(9:17))*10^str2double( file_a(19:21) );  %����ʼʱ��
file_a = fgetl(fid);depth = str2double( file_a(15:18) );range1 = str2double( file_a(24:28) );   %��ȵ����;������
file_a = fgetl(fid);
dz = str2double( file_a(8:16) );   %��������   dz
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
fclose(fid);
%%
range = Depth(:,1);depth = Depth(:,2);
depth_use = interp1( range ,depth , dr:dr:rmax , 'linear');
iz = floor(depth_use/dz);
lengthx = rmax/dr;lengthy = length(aa)/lengthx;
TL_matrc = reshape(aa,lengthy,lengthx);
PR_matrc = reshape(bb,lengthy,lengthx);
PI_matrc = reshape(cc,lengthy,lengthx);
clear aa bb cc;
for i = 1:lengthx
    for j = 1:5
   fjm(i,j) = TL_matrc(iz(i)-3+j,i);
    end
   fjm0(i)  = find( fjm(i) == min( fjm(i) ) );
   fjmpr(i) = PR_matrc(iz(i)-3+fjm0(i),i);
   fjmpi(i) = PI_matrc(iz(i)-3+fjm0(i),i);
   p(i)     = fjmpr(i)+1i*fjmpi(i);
end






