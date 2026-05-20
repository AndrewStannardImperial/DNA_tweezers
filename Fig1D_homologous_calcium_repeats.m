% Andrew Stannard 
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Ca = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];

FRET_mean(1,:) = [0.36 0.96 1.88 4.73 15.17 35.21 51.35 61.56];
FRET_se(1,:) = [0.97 1.58 3.12 3.24 5.74 1.68 0.69 0.84];
FRET_mean(2,:) = [-0.08 1.53 3.48 5.96 17.67 42.62 59.53 69.37];
FRET_se(2,:) = [0.91 0.83 2.00 2.50 5.29 1.91 0.81 0.71];
FRET_mean(3,:) = [-0.64 0.02 0.66 5.11 14.22 42.14 57.66 65.77];
FRET_se(3,:) = [1.17 0.36 1.67 4.11 6.17 0.90 0.95 0.50];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Ca, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(Ca, FRET_mean(1,:), 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(Ca, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(Ca, FRET_mean(2,:), 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(Ca, FRET_mean(3,:), FRET_se(3,:), 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Ca, FRET_mean(3,:), 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)'); xlabel('[Ca^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Ca)/2 max(Ca)*2 -5 75]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-0.12 0.84 2.01 5.27 15.68 39.99 56.18 65.56];
mean_FRET_se = [0.35 0.54 1.00 0.45 1.26 2.93 3.03 2.76];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Ca, mean_FRET_mean, mean_FRET_se, 'Color', CBP(6, :), ...
    'LineWidth', 2);
plot(Ca, mean_FRET_mean, 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(6, :));
ylabel('FRET efficiency (%)'); xlabel('[Ca^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Ca)/2 max(Ca)*2 -5 75]);
legend('', 'average', 'location', 'northwest');