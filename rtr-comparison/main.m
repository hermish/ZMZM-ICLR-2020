%% Parameters
theta = 0.3;
n = 50;
m = 20000;

%% Algorithm
TRIALS = 1;
fileID = fopen('logs50.txt','w');

% CLEAN
times = zeros(TRIALS, 1);
errors = zeros(TRIALS, 1);

for t = 1:TRIALS
    Q = randU(n);
    X = randn(n, m).*(rand(n, m) <= theta);
    Y = Q*X;
    
    tic
    [Q_hat, ~] = DL(Y);
    times(t) = toc;
    errors(t) = 1 - ellFourNorm(Q_hat' * Q) / n;
end

fprintf(fileID, 'Clean (n = %d, p = %d)\n', [n m]);
fprintf(fileID, 'Error: %.2f\\%%\n', mean(errors) * 100);
fprintf(fileID, 'Time: %f\n', mean(times));

% NOISY: 0.2
times = zeros(TRIALS, 1);
errors = zeros(TRIALS, 1);

for t = 1:TRIALS
    Q = randU(n);
    X = randn(n, m).*(rand(n, m) <= theta);
    Y = Q*X + 0.2 * randn(n, m);
    
    tic
    [Q_hat, ~] = DL(Y);
    times(t) = toc;
    errors(t) = 1 - ellFourNorm(Q_hat' * Q) / n;
end

fprintf(fileID, 'Noisy 0.2 (n = %d, p = %d)\n', [n m]);
fprintf(fileID, 'Error: %.2f\\%%\n', mean(errors) * 100);
fprintf(fileID, 'Time: %f\n', mean(times));

% NOISY: 0.4
times = zeros(TRIALS, 1);
errors = zeros(TRIALS, 1);

for t = 1:TRIALS
    Q = randU(n);
    X = randn(n, m).*(rand(n, m) <= theta);
    Y = Q*X + 0.4 * randn(n, m);
    
    tic
    [Q_hat, ~] = DL(Y);
    times(t) = toc;
    errors(t) = 1 - ellFourNorm(Q_hat' * Q) / n;
end

fprintf(fileID, 'Noisy 0.4 (n = %d, p = %d)\n', [n m]);
fprintf(fileID, 'Error: %.2f\\%%\n', mean(errors) * 100);
fprintf(fileID, 'Time: %f\n', mean(times));

% OUTLIERS: 0.2
times = zeros(TRIALS, 1);
errors = zeros(TRIALS, 1);

for t = 1:TRIALS
    Q = randU(n);
    X = randn(n, m).*(rand(n, m) <= theta);
    outliers = randn(n, round(0.2 * m));
    Y = cat(2, Q*X, outliers);
    
    tic
    [Q_hat, ~] = DL(Y);
    times(t) = toc;
    errors(t) = 1 - ellFourNorm(Q_hat' * Q) / n;
end

fprintf(fileID, 'Outliers 0.2 (n = %d, p = %d)\n', [n m]);
fprintf(fileID, 'Error: %.2f\\%%\n', mean(errors) * 100);
fprintf(fileID, 'Time: %f\n', mean(times));

% OUTLIERS: 0.4
times = zeros(TRIALS, 1);
errors = zeros(TRIALS, 1);

for t = 1:TRIALS
    Q = randU(n);
    X = randn(n, m).*(rand(n, m) <= theta);
    outliers = randn(n, round(0.4 * m));
    Y = cat(2, Q*X, outliers);
    
    tic
    [Q_hat, ~] = DL(Y);
    times(t) = toc;
    errors(t) = 1 - ellFourNorm(Q_hat' * Q) / n;
end

fprintf(fileID, 'Outliers 0.4 (n = %d, p = %d)\n', [n m]);
fprintf(fileID, 'Error: %.2f\\%%\n', mean(errors) * 100);
fprintf(fileID, 'Time: %f\n', mean(times));

% CORRUPTIONS: 0.2
times = zeros(TRIALS, 1);
errors = zeros(TRIALS, 1);

for t = 1:TRIALS
    Q = randU(n);
    X = randn(n, m).*(rand(n, m) <= theta);
    Y = Q*X;
    rademacher = (rand(n, m) < 0.5) * 2 - 1;
    Y = Y + (rand(n, m) < 0.2) .* rademacher;
    
    tic
    [Q_hat, ~] = DL(Y);
    times(t) = toc;
    errors(t) = 1 - ellFourNorm(Q_hat' * Q) / n;
end

fprintf(fileID, 'Corruptions 0.2 (n = %d, p = %d)\n', [n m]);
fprintf(fileID, 'Error: %.2f\\%%\n', mean(errors) * 100);
fprintf(fileID, 'Time: %f\n', mean(times));

% CORRUPTIONS: 0.4
times = zeros(TRIALS, 1);
errors = zeros(TRIALS, 1);

for t = 1:TRIALS
    Q = randU(n);
    X = randn(n, m).*(rand(n, m) <= theta);
    Y = Q*X;
    rademacher = (rand(n, m) < 0.5) * 2 - 1;
    Y = Y + (rand(n, m) < 0.4) .* rademacher;
    
    tic
    [Q_hat, ~] = DL(Y);
    times(t) = toc;
    errors(t) = 1 - ellFourNorm(Q_hat' * Q) / n;
end

fprintf(fileID, 'Corruptions 0.4 (n = %d, p = %d)\n', [n m]);
fprintf(fileID, 'Error: %.2f\\%%\n', mean(errors) * 100);
fprintf(fileID, 'Time: %f\n', mean(times));


%% Functions
function value = ellFourNorm(X)
    flattened = reshape(X, [], 1);
    value = norm(flattened, 4) ^ 4;
end