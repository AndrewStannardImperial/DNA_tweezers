% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

conc = [1e-9 2e-9 5e-9 10e-9 20e-9 50e-9 100e-9 200e-9];

FRET_mean(1,:) = [26.94 -1.30 19.12 21.73 18.29 6.85 -0.97 2.77];
FRET_se(1,:) = [2.00 1.08 0.71 0.70 0.54 0.81 1.77 1.18];
FRET_mean(2,:) = [12.08 18.16 16.54 17.41 13.85 8.34 0 0];
FRET_se(2,:) = [2.30 2.50 1.09 0.76 1.40 0.73 0 0];

subplot(121)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(conc, FRET_mean(1,:), FRET_se(1,:), 'Color', CBP(1, :), ...
    'LineWidth', 2);
plot(conc, FRET_mean(1,:), 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(conc, FRET_mean(2,:), FRET_se(2,:), 'Color', CBP(2, :), ...
    'LineWidth', 2);
plot(conc, FRET_mean(2,:), 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
ylabel('FRET efficiency (%)');
xlabel('DNA construct concentration (nM)'); 
xticks([1e-9 10e-9 100e-9]); xticklabels({'1 nM', '10 nM', '100 nM'});
axis([min(conc)/2 max(conc)*2 -5 75]);
legend('', 'repeat 1', '', 'repeat 2', '', 'repeat 3', ...
    'location', 'southwest');

mean_FRET_mean = [19.51 18.16 17.83 19.57 16.07 7.60 1.05 1.31];
mean_FRET_se = [7.43 5.00 1.82 3.05 3.14 1.05 1.09 1.18];

subplot(122)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(conc, mean_FRET_mean, mean_FRET_se, 'Color', CBP(7, :), ...
    'LineWidth', 2);
plot(conc, mean_FRET_mean, 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
ylabel('FRET efficiency (%)');
xlabel('DNA construct concentration (nM)'); 
xticks([1e-9 10e-9 100e-9]); xticklabels({'1 mM', '10 mM', '100 mM'});
axis([min(conc)/2 max(conc)*2 -5 75]);
legend('', 'average', 'location', 'southwest');