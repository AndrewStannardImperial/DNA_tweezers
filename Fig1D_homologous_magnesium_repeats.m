% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Mg = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];

FRET_mean(1,:) = [-1.79 0.22 0.13 2.18 8.07 23.39 40.66 53.11];
FRET_se(1,:) = [0.41 0.80 0.53 1.16 0.81 1.30 0.50 0.69];
FRET_mean(2,:) = [1.13 1.66 2.33 3.75 9.66 24.93 39.62 51.00];
FRET_se(2,:) = [0.78 1.21 0.50 0.45 0.22 0.67 0.42 0.42];
FRET_mean(3,:) = [-0.20 1.76 1.12 2.65 7.40 22.78 38.27 49.18];
FRET_se(3,:) = [0.58 0.56 0.37 0.65 0.30 0.48 1.49 1.18];

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
axis([min(Mg)/2 max(Mg)*2 -5 60]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-0.29 1.21 1.19 2.86 8.38 23.70 39.52 51.10];
mean_FRET_se = [1.03 0.61 0.78 0.57 0.82 0.78 0.85 1.39];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, mean_FRET_mean, mean_FRET_se, 'Color', CBP(7, :), ...
    'LineWidth', 2);
plot(Mg, mean_FRET_mean, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -5 60]);
legend('', 'average', 'location', 'northwest');