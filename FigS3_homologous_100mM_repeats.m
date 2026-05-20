% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

conc = [1e-9 2e-9 5e-9 10e-9 20e-9 50e-9 100e-9 200e-9];

FRET_mean(1,:) = [33.39 38.67 42.09 42.82 39.84 43.62 23.10 13.94];
FRET_se(1,:) = [1.42 0.99 1.17 0.77 1.29 0.46 0.41 1.02];
FRET_mean(2,:) = [42.77 38.07 44.42 41.78 39.48 40.85 17.61 12.83];
FRET_se(2,:) = [0.72 1.25 1.10 0.37 0.49 1.16 0.72 0.74];
FRET_mean(3,:) = [36.05 38.48 39.90 38.82 40.50 41.32 20.14 12.69];
FRET_se(3,:) = [0.93 1.28 5.23 1.14 0.52 0.83 0.98 0.96];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(conc, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(conc, FRET_mean(1,:), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(conc, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(conc, FRET_mean(2,:), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(conc, FRET_mean(3,:), FRET_se(3,:), 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(conc, FRET_mean(3,:), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)');
xlabel('DNA construct concentration (nM)'); 
xticks([1e-9 10e-9 100e-9]); xticklabels({'1 nM', '10 nM', '100 nM'});
axis([min(conc)/2 max(conc)*2 -5 75]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'southwest');

mean_FRET_mean = [37.40 38.41 42.14 41.14 39.94 41.93 20.28 13.15];
mean_FRET_se = [3.42 0.22 1.60 1.47 0.37 1.05 1.94 0.48];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(conc, mean_FRET_mean, mean_FRET_se, 'Color', CBP(7, :), ...
    'LineWidth', 2);
plot(conc, mean_FRET_mean, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
ylabel('FRET efficiency (%)');
xlabel('DNA construct concentration (nM)'); 
xticks([1e-9 10e-9 100e-9]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(conc)/2 max(conc)*2 -5 75]);
legend('', 'average', 'location', 'southwest');