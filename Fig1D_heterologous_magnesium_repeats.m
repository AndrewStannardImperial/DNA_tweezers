% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Mg = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];

FRET_mean(1,:) = [-0.69 0.74 -0.17 1.07 6.48 17.46 29.23 40.63];
FRET_se(1,:) = [0.51 0.64 0.83 0.38 0.43 0.36 0.58 0.50];
FRET_mean(2,:) = [0.79 2.59 2.66 1.89 5.98 16.89 30.28 40.72];
FRET_se(2,:) = [1.14 0.46 0.82 0.10 0.84 1.05 0.94 0.78];
FRET_mean(3,:) = [-3.18 -2.80 -2.43 -0.95 1.30 8.34 20.42 30.40];
FRET_se(3,:) = [0.46 0.64 0.57 0.53 0.41 0.65 0.60 0.60];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(1,:), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(Mg, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(2,:), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(Mg, FRET_mean(3,:), FRET_se(3,:), 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(3,:), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -5 50]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-1.03 0.18 0.02 0.67 4.59 14.23 26.64 37.25];
mean_FRET_se = [1.16 1.58 1.47 0.84 1.65 2.95 3.13 3.43];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, mean_FRET_mean, mean_FRET_se, 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Mg, mean_FRET_mean, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -5 50]);
legend('', 'average', 'location', 'northwest');