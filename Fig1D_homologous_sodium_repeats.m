% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Na = [100e-3 200e-3 500e-3 1 2];

FRET_mean(1,:) = [0.14 -0.05 1.27 2.85 8.99];
FRET_se(1,:) = [0.64 1.46 1.05 0.50 1.32];
FRET_mean(2,:) = [-0.56 0.99 -0.03 3.52 6.41];
FRET_se(2,:) = [0.28 0.52 0.34 0.77 0.81];
FRET_mean(3,:) = [-0.07 0.52 0.70 2.86 5.93];
FRET_se(3,:) = [0.42 0.78 1.17 0.55 0.46];

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
axis([min(Na)/2 max(Na)*2 -2 10]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'northwest');

mean_FRET_mean = [-0.16 0.49 0.65 3.08 7.11];
mean_FRET_se = [0.25 0.37 0.46 0.27 1.16];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Na, mean_FRET_mean, mean_FRET_se, 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(Na, mean_FRET_mean, '<', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
ylabel('FRET efficiency (%)'); xlabel('[Na^+]');
xticks([100e-3 1]); xticklabels({'100 mM', '1 M'});
axis([min(Na)/2 max(Na)*2 -2 10]);
legend('', 'average', 'location', 'northwest');