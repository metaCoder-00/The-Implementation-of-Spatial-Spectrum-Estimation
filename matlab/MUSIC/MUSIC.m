%----------SNR: signal noise radio, n: the number of snapshots, M: the number of array.-------%
SNR = 20;
n = 200;
M = 8;
%----------Consider a ULA, where the array sapcing is a half wavelength of the signal.--------%
c = 3e8;
f = 2.4e9;
lambda = c/f;
d = lambda/2;

%---------Sample: the number of snapshots = n = 200-------------------------------------------%
fs = 3*f;
Ts = 1/fs;
Ns = (0: n - 1)*Ts;

%----------Consider noises are generated from a zero mean Gaussian distribution.--------------%
sigma_N = 0.1;
N = (sigma_N*randn(M, n) + 0).*exp(-1j*2*pi*f*Ns);
%----------Consider three uncorrlated sources at -10 degree, 0 degree and 10 degree,----------------------%
%----------Each source is generated from a zerom mean Gaussian distribution.------------------%
D = 3;                                                         % The number of sources    
sigma_S = sigma_N * 10^(SNR/10);
S = (sigma_S*randn(D, n) + 0).*exp(-1j*2*pi*f*Ns);

theta_S = [-10; 0; 10];
x_k = (0: M - 1)'*d;

A = zeros(M, length(theta_S));
for iter = 1: length(theta_S)
    A(:, iter) = exp(-1j*2*pi*f*(x_k*sind(theta_S(iter)) / c)); % Manifold matrix
end

%---------X = AS + N--------------------------------------------------------------------------%
X = A*S + N;

%---------MUSIC-------------------------------------------------------------------------------%
R_h = (X*X')/n;                                       % R hat: covariance matrix of samples
[V,D] = eig(R_h);                                     % Eigen factorization
D = diag(D);
[D, pin] = sort(D, 'descend');
numOfSrc = AIC(n, M, D);                              % The number of sources
Us = V(:, pin(1: numOfSrc));                          % The subspace of sources
theta = (-90: 0.1: 90)';
Pmu = zeros(length(theta), 1);
for iter = 1: length(theta)
    a = exp(-1j*2*pi*f*(x_k*sind(theta(iter))/c));    % Steering vector
    Pmu(iter) = 1/(a'*(eye(size(Us, 1)) - Us*Us')*a);
end
semilogy(theta, abs(Pmu))
xlabel('angle(degree)')
ylabel('P(dB)')
title(sprintf('SNR = %d, n = %d, M = %d', SNR, n, M))
