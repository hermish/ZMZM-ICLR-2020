function Q = gs(A)
%%% 
% Take a matrix A of full column rank, and returns the result of Gram-Schmidt
% process on the columns of A. 
%%% 

[m, n] = size(A); 
Q = zeros(m, n); 

Q(:, 1) = A(:, 1)/norm(A(:, 1)); 

for ii = 2:n 
	v = A(:, ii); 
	v = v - Q(:, 1:ii-1) * (Q(:, 1:ii-1)' * v);
	Q(:, ii) = v/norm(v);  
end 

