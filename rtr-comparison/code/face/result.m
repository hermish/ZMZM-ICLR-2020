figure;
Compo = 10;
for k = 1:Compo
    x = Y(:,k);
    x = x/max(x);
    x = reshape(x,[192,168]);
    subplot(1,Compo,k);imshow(x);
end
% clc;close all;clear all;
% load data.mat;
x = Y(:,2);
 
%  z = min(x);
%  x = x-z;
%  x = x./max(x);
% x = reshape(x,[192,168]);
% figure; imshow(x);