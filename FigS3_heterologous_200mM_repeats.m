% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

conc = [1e-9 2e-9 5e-9 10e-9 20e-9 50e-9 100e-9];

FRET_mean(1,:) = [31.01 42.86 41.76 44.74 40.55 28.80 20.56];
FRET_se(1,:) = [0.66 0.59 0.47 0.32 0.32 0.82 0.74];
FRET_mean(2,:) = [0 29.78 31.01 30.85 31.79 35.34 0];
FRET_se(2,:) = [0 1.48 0.49 2.01 1.42 1.30 0];
FRET_mean(3,:) = [41.34 37.24 34.60 37.79 40.28 36.83 25.47];
FRET_se(3,:) = [2.17 0.36 0.87 0.36 0.35 0.74 1.06];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(conc, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(conc, FRET_mean(1,:), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(conc(2:6), FRET_mean(2,2:6), FRET_se(2,2:6), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(conc(2:6), FRET_mean(2,2:6), '^', 'MarkerEdgeColor', [0 0 0], ...
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

mean_FRET_mean = [36.18 36.63 35.79 37.79 37.54 33.66 23.02];
mean_FRET_se = [5.16 3.79 3.16 4.01 2.88 2.47 2.46];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(conc, mean_FRET_mean, mean_FRET_se, 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(conc, mean_FRET_mean, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
ylabel('FRET efficiency (%)');
xlabel('DNA construct concentration (nM)'); 
xticks([1e-9 10e-9 100e-9]); xticklabels({'1 nM', '10 nM', '100 nM'});
axis([min(conc)/2 max(conc)*2 -5 75]);
legend('', 'average', 'location', 'southwest');