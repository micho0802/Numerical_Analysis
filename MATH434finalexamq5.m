%Question 5

%a)
%https://www.funnymeetingbackgrounds.com
%/articles/30-fifa-world-cup-qatar-2022-zoom-teams-backgrounds.html

%b)
%covert image to gray scale
clear all, clc
A = imread("MATH434finalexam_q5_data.png");
A_gray = rgb2gray(A);

%c)
%MATLAB code to remove an object from the image

% Convert RGB image to chosen color space
I = rgb2ycbcr(A);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 13.000;
channel1Max = 120.000;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 125.000;
channel2Max = 136.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 147.000;
channel3Max = 182.000;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
remove_objects = A;

% Set background pixels where BW is false to zero.
remove_objects(repmat(~BW,[1 1 3])) = 0;

%d)
%MATLAB code to fill the space of the eliminated objects

fill_back_ground = imfill(BW);

%e)
%MATLAB code to change black table to white and white walls to black

change_bg_wb = imfill(~BW);
imshow(change_bg_wb);