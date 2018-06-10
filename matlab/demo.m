SNR = 25;
snapshots = 200;
sensorNum = 5;
[theta, P] = l1_svd(SNR, snapshots, sensorNum);
subplot(2, 2, 1)
semilogy(theta, abs(P)/max(abs(P)), 'k-', 'LineWidth', 4.0)
set(gca,'XTick',[-100: 25: 100])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10.5)
grid on
xlabel('angle(dgree)')
ylabel('P(dB)')
title(strcat('SNR: ', num2str(SNR), ' snapshots: ', num2str(snapshots), ' sensorNum: ', num2str(sensorNum)))

SNR = 25;
snapshots = 200;
sensorNum = 8;
[theta, P] = l1_svd(SNR, snapshots, sensorNum);
subplot(2, 2, 2)
semilogy(theta, abs(P)/max(abs(P)), 'k-', 'LineWidth', 4.0)
set(gca,'XTick',[-100: 25: 100])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10.5)
grid on
xlabel('angle(dgree)')
ylabel('P(dB)')
title(strcat('SNR: ', num2str(SNR), ' snapshots: ', num2str(snapshots), ' sensorNum: ', num2str(sensorNum)))

SNR = 25;
snapshots = 200;
sensorNum = 16;
[theta, P] = l1_svd(SNR, snapshots, sensorNum);
subplot(2, 2, 3)
semilogy(theta, abs(P)/max(abs(P)), 'k-', 'LineWidth', 4.0)
set(gca,'XTick',[-100: 25: 100])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10.5)
grid on
xlabel('angle(dgree)')
ylabel('P(dB)')
title(strcat('SNR: ', num2str(SNR), ' snapshots: ', num2str(snapshots), ' sensorNum: ', num2str(sensorNum)))

SNR = 25;
snapshots = 200;
sensorNum = 24;
[theta, P] = l1_svd(SNR, snapshots, sensorNum);
subplot(2, 2, 4)
semilogy(theta, abs(P)/max(abs(P)), 'k-', 'LineWidth', 4.0)
set(gca,'XTick',[-100: 25: 100])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10.5)
grid on
xlabel('angle(dgree)')
ylabel('P(dB)')
title(strcat('SNR: ', num2str(SNR), ' snapshots: ', num2str(snapshots), ' sensorNum: ', num2str(sensorNum)))

set(gcf, 'unit', 'inches', 'position', [0.5, 0.5, 10, 10])