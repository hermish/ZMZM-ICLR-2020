clc;close all;clear all;

addpath(genpath('manopt_cur'));
addpath(genpath('test_slate'));
addpath(genpath('image_experiment_focm'));

Times = 100;
patch_size = 8;
D = dir(fullfile('test_slate', '*.pgm'));
idx = 10; 


img = imread(D(idx).name); 

if ndims(img) > 2,
    I = double(rgb2gray(imread(img)));
else
    I = double(img);
end
% figure(7);
% imagesc(img);
Y = image_to_patches(I, patch_size);
[U, S, V] = svd(Y, 'econ'); 
Y = U * V'; 

all_obj = [];

rng(2016);
% Y = (orth(Y'))';
for k = 1:Times
	
	disp(['k = ' num2str(k)]); 
	
	tic 
	
    [A,X] = DL_1(Y);
    all_obj = [all_obj, norm1(inv(A)*Y)];
    
    if k == 1,
        figure(1);
        imagesc(I);
        colormap gray;
        axis off;
        axis image;
        
%        fig = figure(2);
%        set(fig, 'Position', [300, 300, 512, 512]);
%        visualize_orthobasis(A);
    end
    
    figure(3);
    stem(all_obj);
    ylim([0, 1.25 * max(all_obj)]);
    xlabel('Repetition Index'); 
  	ylabel('$$\|\widehat{\mathbf A}^{-1} \mathbf{Y}\|_1$$', 'Interpreter','latex'); 
    
    
    pause(0.1);
    
    toc 
end

    figure(3);
    set(gcf,'PaperPositionMode','auto'); 
    print(gcf, '-dpdf', fullfile('results_new', [num2str(idx), '_obj.pdf'])); 
    system(['pdfcrop ', fullfile('results_new', [num2str(idx), '_obj.pdf'])]); 
