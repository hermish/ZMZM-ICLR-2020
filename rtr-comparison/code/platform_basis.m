clc;close all;clear all;

addpath(genpath('manopt_cur'));
addpath(genpath('test_slate'));
addpath(genpath('image_experiment_focm'));

patch_size = 8;
Basis = cell(3,1);
Image = cell(3,1);
D = dir(fullfile('test_slate', '*.pgm'));
indx = [7,8,10];
for k = 1:3
    img = imread(D(indx(k)).name);
    Image{k} = img;
%     imshow(img);
    if ndims(img) > 2,
        I = double(rgb2gray(imread(img)));
    else
        I = double(img);
    end

    Y = image_to_patches(I, patch_size);
    
    rng(2016);
%     Y = 1/(Y*Y') *Y;
%     Y = (orth(Y'))';
    [U,S,V] = svd(Y);
    S = diag(S).^(-1);
    Pre_cond = U*diag(S)*U';
    Y = Pre_cond*Y;
    
    [A,X] = DL(Y);
    Basis{k} = A;
    
    
    figure;
    imagesc(I);
    colormap gray;
    axis off;
    axis image;
    
    figure;
%     set(fig, 'Position', [300, 300, 512, 512]);
    visualize_orthobasis(A);
    
    
    
    pause(1);
end

