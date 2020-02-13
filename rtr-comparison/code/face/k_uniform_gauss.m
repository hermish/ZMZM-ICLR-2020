function M = k_uniform_gauss( k,m, n )

Omega = zeros(m,n);

for i = 1:n,
    p = randperm(m);
    Omega(p(1:k),i) = 1;
end

V = randn(m,n);

M = Omega .* V;
