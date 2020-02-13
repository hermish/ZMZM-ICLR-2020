% A BFGS solver implemented based on the GRANSO package 1.6:
%   http://www.timmitchell.com/software/GRANSO/

% Update: 
% 10/26/2018 Initial release 

close all; clear all; clc
addpath(genpath('GRANSO'));

%rng(12);
theta = .3;   % sparsity level
n = 300;   % dimension
m = round(10*n^2);    % number of measurements

X = randn(n,m) .* (rand(n,m) <= theta);   % Bernoulli-Gaussian model

% dictionary A=I, one vector (row) q at a time

% supply objective and gradient
obj = @(q) deal(1/m*norm(q'*X, 1), 1/m*X*sign(X'*q));
% supply constraint and gradient
ec = @(q) deal(q'*q - 1, 2*q);

opts.maxit = 2000;
opts.opt_tol = 1e-6;
%opts.limited_mem_size = 250;   % lbfgs
opts.fvalquit = 1e-6;
opts.print_level = 1;
opts.print_frequency = 10; 

tic;
soln = granso(n, obj, [], ec, opts);   % randomly initialize
toc


max(abs(soln.final.x))   % should be close to 1
%stem(soln.final.x)
