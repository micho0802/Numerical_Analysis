clear all, close all, clc
A = imread('roblox4.jpg');
B = rgb2gray(A);
nx = size(B,1); ny = size(B,2);
%% FFT Compression
%Add noise
noise = B + randn(size(B));
Bnoise = fft2(noise);
%Take noise out
Btshiftfilt = Btshift.*ind;
Ffilt = log(abs(Btshiftfilt)+1);
Btfilt = ifftshift(Btshiftfilt); %Apply inverse fft to reomve noise
Bfilt = ifft2(Btfilt);
Bt=fft2(B);    % B is grayscale image from above
Btsort = sort(abs(Bt(:)));  % Sort by magnitude
figure, subplot(2,2,1)
imagesc(B), axis off, colormap gray 
title('Original')
plotind = 2;
% Zero out all small coefficients and inverse transform
for keep=[.001 .002 .003];
    thresh = Btsort(floor((1-keep)*length(Btsort)));
    ind = abs(Bt)>thresh;      % Find small indices
    Atlow = Bt.*ind;           % Threshold small indices
    Alow=uint8(ifft2(Atlow));  % Compressed image
    subplot(2,2,plotind), plotind = plotind + 1;
    imagesc(Alow), axis off 
    ke=keep*100;

    title([num2str(ke,'%2.1f'),'% FFT storage']);
end
set(gcf,'Position',[100 100 550 400])

figure
imagesc(B), axis off, colormap gray 
Image = getframe(gcf);
imwrite(Image.cdata,'100percentFFT.jpg')
for keep=[.001 .002 .003];
thresh = Btsort(floor((1-keep)*length(Btsort)));
    ind = abs(Bt)>thresh;      % Find small indices
    Atlow = Bt.*ind;           % Threshold small indices
    Alow=uint8(ifft2(Atlow));  % Compressed image
imagesc(Alow), axis off 
Image = getframe(gcf);
ke=keep*100;
imwrite(Image.cdata,[num2str(ke,'%2.1f'),'percentFFT.jpg'])
end
close

