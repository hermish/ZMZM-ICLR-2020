clc;close all;clear all;
% Implementation for the Phase Transition of the Trust Region (TR) Method presented in the paper
% "Complete Dictionary Recovery over the Sphere",
% by Ju Sun, Qing Qu, and John Wright.
%
% 1. l1_exp_approx.m: evaluate the value, gradient and hessian
% for the function: h_mu(x) = mu * log cosh(x/mu)
%
% 2. TR_sphere.m: trust region method over the sphere for recovering a single
% sparse vector
%
% 3. Solve_TR_Subproblem.m: solving the trust region subproblem in
% TR_sphere.m using cvx toolbox, MUST BE INSTALLED ahead!
% It can be downloaded from http://cvxr.com/cvx/download/.
%
% Problem formulation: Y = A_0*X_0,
% A_0: n-by-n orthogonal matrix, fixed to  be A_0 = I;
% X_0: n-by-p matrix, each column with fixed k support and standard Gaussian
% distribution for nonzero entries
% We recover one of the sparsest rows in Y = X_0 by solving
%  min_q  h_mu(q' * Y),  s.t.  ||q|| =1.
%
%
% Code written by Ju Sun, Qing Qu and John Wright.
% Last Updated: Sat 25 Apr 2015 11:01:19 PM EDT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment Settings
addpath(genpath('manopt_cur'));
n = 5:10:185; % subspace dimension
p = 10:75:2500;
theta = 0.2;
mu = 1/100;% smoothing parameter mu

num_rep = 5; % number of simulation repetitions
tol = mu; % error tolerance
Prob = zeros(length(p),length(n));% record the recovery probability

rng(1,'twister'); % fix the seed for random number generation
%% Main Simulation
for t = 1:num_rep
    for i = 1:length(p)
        for j = 1:length(n)
            if n(j) > p(i)
                continue;
            end
            % generate the data
            k = ceil( theta*n(j) );
            Y = zeros(n(j),p(i));
            for ll = 1:p(i)
                % generate p columns of k-sparse vectors for Y
                Y(randperm(n(j),k),ll) = randn(k, 1);
            end
            
            % solve the TR-problem for finding one sparse vector
            z_init = randn(n(j),1);
            z_init = z_init/norm(z_init);
            manifold = spherefactory(n(j));
            problem.M = manifold;
            
%             problem.cost = @(q) sum(mu*log(cosh(q'*Y/mu)) );
            problem.cost = @(q)   sum(mu*log( cosh(  (q'*Y/mu).*(abs(q'*Y/mu)<50) ) )....
                +( mu*abs( (q'*Y/mu).*(abs(q'*Y/mu)>50) )  - mu * log(2) ) )  ;
            problem.egrad = @(q)  Y * tanh(Y'*q/mu) ;
            problem.ehess = @(q,u) 1/mu * Y *  ( ( 1- (tanh(Y'*q/mu )).^2 ) .* (Y'*u) ) ;
            
            
            % TRM parameters
            % maxiter can be tuned to be much smaller values when measurements are enough
            options.maxiter = 5000;
            options.tolgradnorm = 1e-4;
            options.verbosity = 0;
            [q,qcost,info,options] = trustregions(problem,z_init, options);
            
            % solve the TR-problem for finding one sparse vector
            
            [~,idx] = max(abs(q));
            ei = zeros(n(j),1);
            ei(idx) = sign(q(idx));
            err = min(norm( q - ei ),norm(q + ei));
            if err<=tol %judging correctness
                Prob(i,j) = Prob(i,j) + 1;
            end
            
            % print intermediate results
            fprintf('Simulation=%d, Sample Number = %d, Dimension = %d, Recovery Error = %f\n'...
                ,t,p(i),n(j),err);
            figure(1);
            imagesc(n, p, Prob/t);
            set(gca,'YDir','normal');
            xlabel('Dimension n');
            ylabel('Sample Number p');
            title('Phase Transition: \theta = 0.2');
            colormap('gray');
            colorbar;
            pause(.25);
        end
    end
    
end

%% Plot Result
figure(1);
imagesc(n, p, Prob/t);
set(gca,'YDir','normal');
xlabel('Dimension n');
ylabel('Sample Number p');
title('Phase Transition: \theta = 0.2');
colormap('gray');
colorbar;
