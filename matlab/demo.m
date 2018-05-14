tic
[theta, P] = l1_svd(25, 200, 8);
% [theta, P] = SpatialSmoothing_MUSIC(25, 300, 8);
semilogy(theta, abs(P)/max(abs(P)))
grid on
[peakVals, locs] = findpeaks(abs(P));
[~, peakIdx] = sort(peakVals, 'descend');
peakIdx = locs(peakIdx(1: 3));
theta_hat = theta(peakIdx)
toc