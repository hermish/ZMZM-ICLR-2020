function M = k_uniform_rademacher(k, m, n)

Omega = zeros(m,n);

for i = 1:n,
    p = randperm(m);
    Omega(p(1:k),i) = 1;
end

rand_mtx = rand(m,n) - 0.5;
V = (rand_mtx > 0) - (rand_mtx < 0); 

M = Omega .* V;
