% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Mg = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];

FRET_mean(1,:) = [-0.37 0.27 0.19 -0.34 2.03 8.14 15.53 21.66];
FRET_se(1,:) = [0.69 1.06 1.03 0.51 0.62 0.99 0.84 1.87];
FRET_mean(2,:) = [-1.04 -1.12 0.67 0.71 1.70 7.66 15.59 21.73];
FRET_se(2,:) = [0.89 1.33 0.49 1.35 1.29 1.53 1.99 0.82];
FRET_mean(3,:) = [-0.41 1.72 0.02 -0.08 2.91 10.43 13.37 18.32];
FRET_se(3,:) = [0.75 1.27 0.62 0.33 0.72 1.50 1.79 4.83];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(1,:), 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(Mg, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(2,:), 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(Mg, FRET_mean(3,:), FRET_se(3,:), 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Mg, FRET_mean(3,:), 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -5 25]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-0.61 0.29 0.29 0.10 2.21 8.74 14.83 20.57];
mean_FRET_se = [0.27 1.00 0.24 0.39 0.44 1.05 0.89 1.38];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, mean_FRET_mean, mean_FRET_se, 'Color', CBP(7, :), ...
    'LineWidth', 2);
plot(Mg, mean_FRET_mean, 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
ylabel('FRET efficiency (%)'); xlabel('[Mg^{2+}]');
xticks([1e-3 10e-3 100e-3]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(Mg)/2 max(Mg)*2 -5 25]);
legend('', 'average', 'location', 'northwest');