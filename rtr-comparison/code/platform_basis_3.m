clc;close all;clear all;

addpath(genpath('manopt_cur'));
addpath(genpath('test_slate'));
addpath(genpath('image_experiment_focm'));


rng(2016); 

patch_size = 8;
Basis = cell(3,1);
Image = cell(3,1);
D = dir(fullfile('test_slate', '*.pgm'));
indx = [7,8,10];
for k = 3:3
    img = imread(D(indx(k)).name);
    Image{k} = img;
%     imshow(img);
    if ndims(img) > 2,
        I = double(rgb2gray(imread(img)));
    else
        I = double(img);
    end

    Y = image_to_patches(I, patch_size);

		[U, S, V] = svd(Y, 'econ'); 
		Y = U * V'; 
    
    [A,X] = DL_1(Y);
    Basis{k} = A;
    
    
    figure;
    imagesc(I);
    colormap gray;
    axis off;
    axis image;
    
		fig = figure(2);
		set(fig, 'Position', [300, 300, 512, 512]); 
		visualize_orthobasis(A);
    set(gcf,'PaperPositionMode','auto'); 
    print(gcf, '-dpdf', fullfile('results_new', [num2str(k), '_dict.pdf'])); 
    system(['pdfcrop ', fullfile('results_new', [num2str(k), '_dict.pdf'])]); 
    
end

