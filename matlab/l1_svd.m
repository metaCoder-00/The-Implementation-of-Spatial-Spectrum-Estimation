function [theta, S_sv_hat] = l1_svd(SNR, snapshots, sensorNum)

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
    theta_S = [-25; 0; 25];
    sourceNum = length(theta_S);
    sigma_S = sigma_N * 10^(SNR/10);
    signalCovMat = [sigma_S, 0.99*sigma_S, 0; 0.99*sigma_S, sigma_S, 0; 0, 0, sigma_S];
    signalAmp = mvnrnd(zeros(sourceNum, 1), signalCovMat, snapshots);
    signalPhase = exp(-1j*2*pi*f*Ns);
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

    [~, ~, rightVec] = svd(arrayOut);
    dimReduceMat = [eye(sourceNum), zeros(sourceNum, snapshots - sourceNum)]';
    arrayOut_DimReduce = arrayOut * rightVec * dimReduceMat;
    theta = (-90: 0.1: 90)';                % angle scan range
    manifoldMat_hat = zeros(sensorNum, length(theta));
    for col = 1: length(theta)
        manifoldMat_hat(:, col) = exp(-1j*2*pi*f*((spacingK*sind(theta(col)))/c));
    end

    %----------SOC programming-----------------------%
    regParam = 1.42;
    sumVec = ones(length(theta), 1);
    cvx_begin
        variables p q r(length(theta));
        variable S_sv(length(theta), sourceNum) complex;
        for row = 1: length(theta)
            S_sv_hat(row, :) = norm(S_sv(row, :));
        end
        minimize(p + regParam*q);
        subject to
             Z_k = arrayOut_DimReduce - manifoldMat_hat*S_sv;
             norm(Z_k, 'fro') <= p;

            sumVec' * r <= q;

            for idx = 1: length(theta)
                S_sv_hat(idx) <= r(idx);
            end
    cvx_end
end





