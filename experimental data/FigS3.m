% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

hetero_conc = [1e-9 2e-9 5e-9 10e-9 20e-9 50e-9 100e-9];
hetero_FRET_mean = [36.18 36.63 35.79 37.79 37.54 33.66 23.02];
hetero_FRET_se = [5.16 3.79 3.16 4.01 2.88 2.47 2.46];
homo_conc = [1e-9 2e-9 5e-9 10e-9 20e-9 50e-9 100e-9 200e-9];
homo_FRET_mean = [37.40 38.41 42.14 41.14 39.94 41.93 20.28 13.15];
homo_FRET_se = [3.42 0.22 1.60 1.47 0.37 1.05 1.94 0.48];
homo_long_conc = [1e-9 2e-9 5e-9 10e-9 20e-9 50e-9 100e-9 200e-9];
homo_long_FRET_mean = [19.51 18.16 17.83 19.57 16.07 7.60 1.05 1.31];
homo_long_FRET_se = [7.43 5.00 1.82 3.05 3.14 1.05 1.09 1.18];

set(gca, 'FontSize', 12, 'XScale', 'log'); hold on; axis square;
errorbar(hetero_conc, hetero_FRET_mean, hetero_FRET_se, 'Color', ...
    CBP(3, :), 'LineWidth', 2);
plot(hetero_conc, hetero_FRET_mean, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
errorbar(homo_conc, homo_FRET_mean, homo_FRET_se, 'Color', CBP(7, :), ...
    'LineWidth', 2);
plot(homo_conc, homo_FRET_mean, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
errorbar(homo_long_conc, homo_long_FRET_mean, homo_long_FRET_se, ...
    'Color', CBP(7, :), 'LineWidth', 2);
plot(homo_long_conc, homo_long_FRET_mean, 'd', 'MarkerEdgeColor', ...
    [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
xlabel('DNA construct concentration (nM)'); 
ylabel('Average FRET efficiency, \langleE\rangle (%)');
xticks([1e-9 10e-9 100e-9]); xticklabels({'1', '10', '100'});
axis([min(homo_conc)/2 max(homo_conc)*2 -5 75]);