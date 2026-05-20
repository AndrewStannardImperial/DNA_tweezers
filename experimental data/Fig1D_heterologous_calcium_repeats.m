% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Ca = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];

FRET_mean(1,:) = [-1.21 0.63 0.32 0.43 5.64 25.26 43.24 59.63];
FRET_se(1,:) = [0.59 0.63 0.80 0.54 0.86 0.46 0.32 0.27];
FRET_mean(2,:) = [1.14 2.95 4.85 4.80 11.81 28.85 45.95 51.61];
FRET_se(2,:) = [0.87 0.68 0.33 0.91 0.57 1.03 1.00 2.10];
FRET_mean(3,:) = [-3.06 -1.49 -3.24 -0.07 10.05 24.08 37.24 45.57];
FRET_se(3,:) = [0.23 0.60 0.68 0.75 1.40 0.48 2.02 0.18];

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
axis([min(Ca)/2 max(Ca)*2 -10 65]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-1.04 0.70 0.64 1.72 9.17 26.06 42.14 52.27];
mean_FRET_se = [1.49 1.57 2.86 1.89 2.25 1.76 3.15 4.99];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Ca, mean_FRET_mean, mean_FRET_se, 'Color', CBP(4, :), ...
    'LineWidth', 2);
plot(Ca, mean_FRET_mean, 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(4, :));
ylabel('FRET efficiency (%)'); xlabel('[Ca^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Ca)/2 max(Ca)*2 -10 65]);
legend('', 'average', 'location', 'northwest');