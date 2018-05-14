theta_S = [-20; 0; 20];
SNR = 25;
sensorNum = 8;

N = 1000;

snapshotsVals = (100: 100: 500)';
SSm_MUSIC_RMES = zeros(length(snapshotsVals), 1);
l1_SVD_RMES = zeros(length(snapshotsVals), 1);

for idx = 1: length(snapshotsVals)
    RMES = 0;   
    for itr = 1: N
        [theta, P] = SpatialSmoothing_MUSIC(SNR, snapshotsVals(idx), sensorNum);
        [peakVals, locs] = findpeaks(abs(P));
        if length(locs) < length(theta_S)
            theta_hat = zeros(length(theta_S), 1);
        else
            [~, peakIdx] = sort(peakVals, 'descend');
            peakIdx = locs(peakIdx(1: length(theta_S)));
            theta_hat = theta(peakIdx);
            theta_hat = sort(theta_hat);
        end
        
        res = theta_S - theta_hat;
        RMES = RMES + sqrt((res'*res)/length(theta_S));
    end
    RMES = RMES / N;
    SSm_MUSIC_RMES(idx) = RMES;
    
    RMES = 0;
    for itr = 1: N
        [theta, P] = l1_svd(SNR, snapshotsVals(idx), sensorNum);
        [peakVals, locs] = findpeaks(abs(P));
        if length(locs) < length(theta_S)
            theta_hat = zeros(length(theta_S), 1);
        else
            [~, peakIdx] = sort(peakVals, 'descend');
            peakIdx = locs(peakIdx(1: length(theta_S)));
            theta_hat = theta(peakIdx);
            theta_hat = sort(theta_hat);
        end
        
        res = theta_S - theta_hat;
        RMES = RMES + sqrt((res'*res)/length(theta_S));
    end
    RMES = RMES / N;
    l1_SVD_RMES(idx) = RMES;   
end

figure(2)

plot(snapshotsVals, SSm_MUSIC_RMES, 'o-', snapshotsVals, l1_SVD_RMES, '*--')
legend('Spatial Smoothing MUSIC', 'l1-SVD');
xlabel('The number of snapshots')
ylabel('RMSE')
title(strcat('The number of test: ', num2str(N)))