clear all, close all, clc

A=imread('roblox4.jpg');
X=double(rgb2gray(A));  % Convert RBG to gray, 256 bit to double.
nx = size(X,1); ny = size(X,2);

noise = A + 200*randn(size(A));

[U,S,V] = svd(X);

figure, subplot(2,2,1)
imagesc(X), axis off, colormap gray 
title('Original')
plotind = 2;
for r=[13 14 15];  % Truncation value
    Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Approx. image
    subplot(2,2,plotind), plotind = plotind + 1;
    imagesc(Xapprox), axis off %Display image with scaled colors
    title(['r=',num2str(r,'%d'),', ',num2str(100*r*(nx+ny)/(nx*ny),'%2.2f'),'% SVD storage']);
end
set(gcf,'Position',[100 100 550 400])

figure, 
imagesc(X), axis off, colormap gray 
Image = getframe(gcf);
imwrite(Image.cdata,'100percentSVD.jpg')
for r=[13 14 15];  % Truncation value
    Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Approx. image
    imagesc(Xapprox), axis off %Display image with scaled colors
    Image = getframe(gcf);
    imwrite(Image.cdata,[num2str(100*r*(nx+ny)/(nx*ny),'%2.2f'),'percentSVD.jpg'])
end
close
