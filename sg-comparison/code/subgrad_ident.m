
% Riemannian subgradient descent (identity dictionary)
% Update: 
% 10/26/2018 Initial release 

clear all;

% generate data
theta = .3;   % sparsity level
n = 70;   % dimension
p = [0.5,1,1.5,2,2.5];   % sample complexity (as power of n)
max_iter = 3500;
eps = 10^-3;   % neighborhood size (stopping criteria)

success_vec = zeros(length(p), 1);

for k = 1:length(p)
	m = round(10*n^p(k));    % number of measurements
	success = zeros(10, 1);

	for j = 1:10  % 10 random instances
		X = randn(n, m).*(rand(n, m) <= theta);   % iid Bern-Gaussian model
		record_basis = zeros(n, 1);

		for l = 1:round(5*n*log(n))
			disp(['k = ', num2str(k), ', j = ', num2str(j), ', l = ', num2str(l)]);

			tic
			% random initialization on sphere
			q = randn(n,1);
			q = q./norm(q);

			for i = 1:max_iter
				eta = 1/sqrt(i);  % diminishing step size
				egrad = 1/m*X*sign(X'*q);
				q = q - eta * (egrad - q*(q'*egrad));    % Riemannian step
				q = q/norm(q);   % projection

				% early stopping if recovered
				if abs(max(abs(q))-1) <= eps
					[~,ind] = max(abs(q));
					record_basis(ind) = 1;  % recovery up to sign
					break
				end

			end
			toc

			if nnz(~(record_basis(:))) == 0   % found all standard basis
				success(j) = 1;
				break
			end

		end  % end for 5*n*log(n) runs

	end  % end for the instances

	disp(['Success rate for n = ', num2str(n), 'm = ', num2str(m), ': ', num2str(mean(success))]);
	success_vec(k) = mean(success);

end  % end for different powers
success_vec
