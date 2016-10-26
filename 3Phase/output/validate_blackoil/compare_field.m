clear all;
close all;
NX = 50;
NY = 50;

load results.dat
load ecl_rs.dat;
load ecl_p.dat;
load ecl_sw.dat;
load ecl_sg.dat;
p = reshape( results(:,1),NX,NY);
sw = reshape( results(:,2),NX,NY);
sg = reshape( results(:,3),NX,NY);
so = reshape( results(:,4),NX,NY);
Rs = reshape( results(:,5),NX,NY);
ecl_p=reshape( ecl_p',NX,NY);
clear results
figure;
subplot(1,2,1),imagesc(ecl_p),title('eclipse pressure');
axis square
colorbar
subplot(1,2,2),imagesc(p),colormap('jet'),title('pressure');
axis square
colorbar
ecl_sw=reshape( ecl_sw',NX,NY);
clear results
figure;
subplot(1,2,1),imagesc(ecl_sw),title('eclipse Sw');
axis square
colorbar
subplot(1,2,2),imagesc(sw),colormap('jet'),title('Sw');
axis square
colorbar
ecl_sg=reshape( ecl_sg',NX,NY);
clear results
figure;
subplot(1,2,1),imagesc(ecl_sg),title('eclipse Sg');
axis square
colorbar
subplot(1,2,2),imagesc(sg),colormap('jet'),title('Sg');
axis square
colorbar
ecl_rs=reshape( ecl_rs',NX,NY);
clear results
figure;
subplot(1,2,1),imagesc(ecl_rs),title('eclipse Rs');
axis square
colorbar
subplot(1,2,2),imagesc(Rs*178.1076),colormap('jet'),title('Rs');
axis square
colorbar