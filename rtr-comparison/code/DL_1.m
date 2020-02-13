function [A,X] = DL(Y)
mu = 1e-2;
[n, p] = size(Y);
N = n;
% [V,D] = eig(Y'*Y);
% S = D^(-1/2);
% Y = (V'*S*V)*Y;
Y_0 = Y;
Q = []; X = [];
U = eye(n);
for t = 1:N
   
    
    z_init = randn(n,1);
    z_init = z_init/norm(z_init);
    manifold = spherefactory(n);
    problem.M = manifold;
    
%     problem.cost = @(q) sum(mu*log(cosh(q'*Y/mu)) );
    problem.cost = @(z)   sum(mu*log( cosh(  (z'*Y/mu).*(abs(z'*Y/mu)<50) ) )....
        +( mu*abs( (z'*Y/mu).*(abs(z'*Y/mu)>50) )  - mu * log(2) ) )  ;
    problem.egrad = @(z)  Y * tanh(Y'*z/mu) ;
    problem.ehess = @(z,u) 1/mu * Y *  ( ( 1- (tanh(Y'*z/mu )).^2 ) .* (Y'*u) ) ;
    options.maxiter = 1e4;
    options.tolgradnorm = 1e-7;
    options.verbosity = 0;
    [z,zcost,info,options] = trustregions(problem,z_init, options);
    r = U * z;
    
%    LP rounding 
    cvx_begin quiet 
    variable q(N, 1);
    minimize (norm(q'*Y_0 , 1));
    subject to
       q'* r == 1;
%     (r'*Y ) * x' == 1;
    cvx_end
    X = [X; q'*Y_0];
    
    Q = [Q, q];
    n = n - 1; 
    U = null(Q'); 
    Y = U' * Y_0;
   
end

A = Y_0 / X;

end

