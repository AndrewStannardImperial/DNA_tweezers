% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Mg = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];

FRET_mean(1,:) = [-0.01 0.95 1.21 1.85 7.75 21.20 41.44 52.91];
FRET_se(1,:) = [1.96 1.59 1.32 1.14 0.87 3.35 0.81 0.96];
FRET_mean(2,:) = [0.50 0.42 2.06 3.14 8.42 24.73 36.54 46.83];
FRET_se(2,:) = [1.25 1.32 1.39 1.15 0.98 0.09 3.41 1.29];
FRET_mean(3,:) = [0.39 -0.32 -0.47 2.11 7.80 19.60 43.86 53.95];
FRET_se(3,:) = [0.46 0.38 0.38 0.94 0.61 3.65 0.38 0.24];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(1,:), 'o', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(Mg, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(2,:), 'o', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(Mg, FRET_mean(3,:), FRET_se(3,:), 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(3,:), 'o', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -5 60]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [0.29 0.35 0.93 2.37 7.99 21.84 40.61 51.23];
mean_FRET_se = [0.19 0.45 0.91 0.48 0.26 1.86 2.64 2.72];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, mean_FRET_mean, mean_FRET_se, 'Color', CBP(5, :), ...
    'LineWidth', 2);
plot(Mg, mean_FRET_mean, 'o', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(5, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -5 60]);
legend('', 'average', 'location', 'northwest');