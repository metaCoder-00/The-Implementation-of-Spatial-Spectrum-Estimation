function [theta, P] = SpatialSmoothing_MUSIC(SNR, snapshots, sensorNum)

    % SNR = 20;
    % snapshots = 200;
    % sensorNum = 8;

    %----------Consider a ULA, where the array sapcing is a half wavelength of the signal.--------%
    c = 3e8;
    f = 2.4e9;                       % frequency is 2.4GHz
    wavelength = c / f;              % lambda
    spacing = wavelength / 2;        % ULA's spacing

    %---------Sample: sample frequency is fs = 3f-------------------------------------------%
    fs = 3 * f;                      
    Ts = 1 / fs;                             % Sample period
    Ns = Ts*(0: snapshots - 1)';          % Sample spacing

    %----------Consider noises are generated from a zero mean Gaussian distribution.--------
    sigma_N = 0.1;
    noiseCovMat = sigma_N * eye(sensorNum);

    noiseAmp = mvnrnd(zeros(sensorNum, 1), noiseCovMat, snapshots);
    noisePhase = exp(-1j*2*pi*f*Ns + randn());
    noiseMat = zeros(size(noiseAmp));               % Each row is A sample 
    for col = 1: sensorNum
        noiseMat(:, col) = noiseAmp(:, col) .* noisePhase;
    end

    %----------Consider three sources at -10 degree, 0 degree and 10 degree.----------------------%
    %----------Sources at -10 degree and 0 degree is coherent.------------------------------------%
    %----------Each source is generated from a zerom mean Gaussian distribution.------------------%
    theta_S = [-10; 0; 10];
    sourceNum = length(theta_S);
    sigma_S = sigma_N * 10^(SNR/20);
    signalCovMat = [sigma_S, 0.99*sigma_S, 0; 0.99*sigma_S, sigma_S, 0; 0, 0, sigma_S];
    signalAmp = mvnrnd(zeros(sourceNum, 1), signalCovMat, snapshots);
    signalPhase = exp(-1j*2*pi*f*Ns + randn());
    signalMat = zeros(size(signalAmp));             % Each row is A sample 
    for col = 1: sourceNum
        signalMat(:, col) = signalAmp(:, col) .* signalPhase;
    end

    spacingK = spacing * (0: sensorNum - 1)';
    manifoldMat = zeros(sensorNum, sourceNum);
    for col = 1: sourceNum
        manifoldMat(:, col) = exp(-1j*2*pi*f*((spacingK*sind(theta_S(col)))/c));
    end
    arrayOut = manifoldMat*signalMat.' + noiseMat.';

    covMat = (arrayOut*arrayOut') / snapshots;
    m = round((sensorNum + (sourceNum/2) + 1)/2);       % m S.T. D <= m <= M - D/2 + 1
    p = sensorNum + 1 - m;
    R_fk = zeros(m, m);
    R_bk = zeros(m, m);

    for idx = 1: p
        Z_k = [zeros(m, idx - 1), eye(m), zeros(m, p - idx)];
        R_fk = R_fk + Z_k*covMat*Z_k';

        Q_k = [zeros(m, idx - 1), fliplr(eye(m)), zeros(m, p - idx)];
        R_bk = R_bk + Q_k*conj(covMat)*Q_k';
    end

    R_f = R_fk / p;
    R_b = R_bk / p;
    R_fb = (R_f + R_b)/2;

    [eigenVec, eigenVals] = eig(R_fb);
    eigenVals = diag(eigenVals);
    [~, eigenValsIdx] = sort(eigenVals);
    noiseSubspace = eigenVec(:, eigenValsIdx(1: m - sourceNum));

    theta = (-90: 90)';
    P = zeros(length(theta), 1);
    for idx = 1: length(theta)
        steeringVec = exp(-1j*2*pi*f*(spacing*(0: m - 1)'*sind(theta(idx)))/c);
        P(idx) = 1 / (steeringVec'*(noiseSubspace*noiseSubspace')*steeringVec);
    end
    
end



