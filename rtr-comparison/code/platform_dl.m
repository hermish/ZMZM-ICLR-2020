clc;close all;clear all;

addpath(genpath('manopt_cur'));
addpath(genpath('test_slate'));
addpath(genpath('image_experiment_focm'));

Times = 100;
patch_size = 8;
D = dir(fullfile('test_slate', '*.pgm'));
img = imread(D(8).name);
% imshow(img);
indx = [7,8,10];
if ndims(img) > 2,
    I = double(rgb2gray(imread(img)));
else
    I = double(img);
end
% figure(7);
% imagesc(img);
Y = image_to_patches(I, patch_size);
all_obj = [];

rng(2016);
% Y = (orth(Y'))';
for k = 1:Times
    [A,X] = DL(Y);
    all_obj = [all_obj, norm1(A'*Y)];
    
    if k == 1,
        figure(1);
        imagesc(I);
        colormap gray;
        axis off;
        axis image;
        
        fig = figure(2);
        set(fig, 'Position', [300, 300, 512, 512]);
        visualize_orthobasis(A);
    end
    
    figure(3);
    stem(all_obj);
    ylim([0, 1.25 * max(all_obj)]);
    
    pause(0.1);
end
