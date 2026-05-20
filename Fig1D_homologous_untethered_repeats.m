% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Mg = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];

FRET_mean(1,:) = [0.15 0.88 1.48 0.91 1.25 1.23 2.33 3.07];
FRET_se(1,:) = [1.01 0.49 0.84 1.07 1.02 1.37 0.71 0.77];
FRET_mean(2,:) = [-1.01 -0.72 -0.73 -1.55 0.90 1.16 1.29 -1.21];
FRET_se(2,:) = [0.81 0.90 0.80 0.31 0.61 0.88 1.18 1.32];
FRET_mean(3,:) = [-0.77 0.00 -0.17 -0.86 -0.22 0.01 0.12 0.54];
FRET_se(3,:) = [0.57 0.78 0.78 0.75 0.88 0.90 0.50 0.63];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(1,:), 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(Mg, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(2,:), 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(Mg, FRET_mean(3,:), FRET_se(3,:), 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(3,:), 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -2 4]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-0.54 0.05 0.19 -0.50 0.64 0.80 1.25 0.80];
mean_FRET_se = [0.43 0.57 0.81 0.90 0.54 0.48 0.78 1.24];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, mean_FRET_mean, mean_FRET_se, 'Color', CBP(7, :), ...
    'LineWidth', 2);
plot(Mg, mean_FRET_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -2 4]);
legend('', 'average', 'location', 'northwest');