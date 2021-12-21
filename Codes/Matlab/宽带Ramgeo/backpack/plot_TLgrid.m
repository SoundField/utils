clear;clc
filename = 'tl.grid';

rd = 0 : 2 : 800;
rr = 10*1 : 10*1 : 52000;

Nrd = length( rd );
Nrr = length( rr );


fid = fopen( filename, 'rb' );

if ( fid == -1 )
   error( 'No shade file with that name exists; you must run a model first' );
end
a  = fread(fid,1,'int64');
TL = fread( fid, [ Nrd, Nrr ], 'float32' );    %Read complex data

fclose( fid );
whos
figure;
pcolor( rr/1e3, rd, TL );
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
shading interp;
jeta=jet(48);  jetb=flipud(jeta);  colormap(jetb);
caxis([30,120]);  
colorbar('NorthOutside')
% ylim([10 120]);
% xlim([0 30]);
