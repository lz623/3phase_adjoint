clear all 
close all
NX = 50;
NY = 50;
% 
% load field.out
% phi =reshape( field(:,1),NX,NY);
% kx = reshape( field(:,2),NX,NY);
% ky = reshape( field(:,3),NX,NY);
% kz = reshape( field(:,4),NX,NY);
% figure;
% imagesc(phi),colormap('jet');
% axis square
% figure;
% imagesc(kx),colormap('jet');
% axis square
% figure;
% imagesc(ky),colormap('jet');
% figure;
% imagesc(kz),colormap('jet');
% axis square
% clear field

load results.dat
p = reshape( results(:,1),NX,NY);
sw = reshape( results(:,2),NX,NY);
sg = reshape( results(:,3),NX,NY);
so = reshape( results(:,4),NX,NY);
Rs = reshape( results(:,5),NX,NY);
clear results
figure;
subplot(1,2,1),imagesc(ecl_p),title('eclipse pressure');
axis square
colorbar
subplot(1,2,1),imagesc(p),colormap('jet'),title('pressure');
axis square
colorbar
figure;
imagesc(sw),colormap('jet'),title('Sw');
axis square
colorbar
figure;
imagesc(sg),colormap('jet'),title('Sg');
axis square
colorbar
figure;
imagesc(so),colormap('jet'),title('So');
axis square
colorbar
figure;
imagesc(Rs),colormap('jet'),title('Ro');
axis square
colorbar
