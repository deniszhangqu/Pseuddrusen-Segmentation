clc;close all; clear all;

% read in the image.
% img = imread('D:\Master-Dataset\RPD3.png');
load('D:\Master-Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1001.mat');
img=images(:,:,50);
%%
% error checking, get one channel from image.
if size(img,3) > 1
    img = rgb2gray(img);
end
% make image type as double.
img = double(img);

% get retinal layers.
[retinalLayers, params] = getRetinalLayers(img);  %函数里面的参数可以设置

