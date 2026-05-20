% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Na = [100e-3 200e-3 500e-3 1 2];

FRET_mean(1,:) = [-1.51 -2.26 -2.65 1.65 7.24];
FRET_se(1,:) = [0.78 0.62 0.87 0.28 0.58];
FRET_mean(2,:) = [0.82 3.07 3.76 6.22 10.06];
FRET_se(2,:) = [0.73 0.45 0.83 0.56 0.59];
FRET_mean(3,:) = [0.45 0.21 0.06 2.25 7.36];
FRET_se(3,:) = [0.61 0.88 2.63 0.95 0.26];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Na, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(Na, FRET_mean(1,:), '<', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(Na, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(Na, FRET_mean(2,:), '<', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(Na, FRET_mean(3,:), FRET_se(3,:), 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Na, FRET_mean(3,:), '<', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)'); xlabel('[Na^+]');
xticks([100e-3 1]); xticklabels({'100 mM', '1 M'});
axis([min(Na)/2 max(Na)*2 -4 12]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-0.08 0.34 0.39 3.37 8.22];
mean_FRET_se = [0.89 1.89 2.28 1.76 1.13];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Na, mean_FRET_mean, mean_FRET_se, 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(Na, mean_FRET_mean, '<', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
ylabel('FRET efficiency (%)'); xlabel('[Na^+]');
xticks([100e-3 1]); xticklabels({'100 mM', '1 M'});
axis([min(Na)/2 max(Na)*2 -4 12]);
legend('', 'average', 'location', 'northwest');