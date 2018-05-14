theta_S = [-20; 0; 20];
SNR = 25;
snapshots = 200;

N = 1000;

sensorNumVals = (6: 4: 34)';
SSm_MUSIC_RMES = zeros(length(sensorNumVals), 1);
l1_SVD_RMES = zeros(length(sensorNumVals), 1);
for idx = 1: length(sensorNumVals)
    RMES = 0;   
    for itr = 1: N
        [theta, P] = SpatialSmoothing_MUSIC(SNR, snapshots, sensorNumVals(idx));
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
        [theta, P] = l1_svd(SNR, snapshots, sensorNumVals(idx));
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


figure(3)

plot(sensorNumVals, SSm_MUSIC_RMES, 'o-', sensorNumVals, l1_SVD_RMES, '*--')
legend('Spatial Smoothing MUSIC', 'l1-SVD');
xlabel('The number of sensor')
ylabel('RMSE')
title(strcat('The number of test: ', num2str(N)))