clear;clc;
unix([pwd,'/rambd']);
filename='tl.grid';
fid = fopen( filename, 'r' );

if ( fid == -1 )
   error( 'No shade file with that name exists; you must run a model first' );
end
r = 10:10:52000;
Nrr = 5200;
% rec=fread(fid,1,'int');
Nff = 32;
TL = fread( fid,[Nrr,Nff],'float' );    %Read complex data
% rec=fread(fid,1,'int');
% TL(2,:) = fread( fid,Nrr,'float' );
fclose(fid);

TL = TL';
%%
% figure;
% for i = 1:Nff
%     plot(r/1e3,TL(i,:));
%     axis ij;hold on;
% end
tl = sum(TL)/Nff;
plot(r/1e3,tl);
axis ij;grid on;