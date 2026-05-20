% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Mg = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];
Ca = [1e-3 2e-3 5e-3 10e-3 20e-3 50e-3 100e-3 200e-3];
Na = [100e-3 200e-3 500e-3 1e0 2e0];

Mg_hetero = [-1.03 0.18 0.02 0.67 4.59 14.23 26.64 37.25];
Mg_hetero_se = [1.16 1.58 1.47 0.84 1.65 2.95 3.13 3.43];
Ca_hetero = [-1.04 0.70 0.64 1.72 9.17 26.06 42.14 52.27];
Ca_hetero_se = [1.49 1.57 2.86 1.89 2.25 1.76 3.15 4.99];
Na_hetero = [-0.08 0.34 0.39 3.37 8.22];
Na_hetero_se = [0.89 1.89 2.28 1.76 1.13];
untethered_hetero = [-0.16 0.28 0.30 0.58 0.23 -0.02 0.46 -0.25];
untethered_hetero_se = [0.19 0.47 0.58 0.27 0.51 0.73 0.85 1.47];
long_hetero = [-0.20 0.68 0.83 0.93 2.52 5.44 10.97 15.20];
long_hetero_se = [0.95 0.35 1.42 0.56 1.19 0.64 1.68 1.34];

Mg_homo = [-0.29 1.21 1.19 2.86 8.38 23.70 39.52 51.10];
Mg_homo_se = [1.03 0.61 0.78 0.57 0.82 0.78 0.85 1.39];
Ca_homo = [-0.12 0.84 2.01 5.27 15.68 39.99 56.18 65.56];
Ca_homo_se = [0.35 0.54 1.00 0.45 1.26 2.93 3.03 2.76];
Na_homo = [-0.16 0.49 0.65 3.08 7.11];
Na_homo_se = [0.25 0.37 0.46 0.27 1.16];
untethered_homo = [-0.54 0.05 0.19 -0.50 0.64 0.80 1.25 0.80];
untethered_homo_se = [0.43 0.57 0.81 0.90 0.54 0.48 0.78 1.24];
long_homo = [-0.61 0.29 0.29 0.10 2.21 8.74 14.83 20.57];
long_homo_se = [0.27 1.00 0.24 0.39 0.44 1.05 0.89 1.38];
control_homo = [0.29 0.35 0.93 2.37 7.99 21.84 40.61 51.23];
control_homo_se = [0.19 0.45 0.91 0.48 0.26 1.86 2.64 2.72];

subplot(121)
set(gca, 'FontSize', 12, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, long_hetero, long_hetero_se, 'Color', CBP(3, :), ...
    'LineWidth', 2);
plot(Mg, long_hetero, 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
errorbar(Mg, untethered_hetero, untethered_hetero_se, 'Color', ...
    CBP(3, :), 'LineWidth', 2);
plot(Mg, untethered_hetero, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
errorbar(Na, Na_hetero, Na_hetero_se, 'Color', CBP(1, :), 'LineWidth', 2);
plot(Na, Na_hetero, '<', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(1, :));
errorbar(Ca, Ca_hetero, Ca_hetero_se, 'Color', CBP(4, :), 'LineWidth', 2);
plot(Ca, Ca_hetero, 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(4, :));
errorbar(Mg, Mg_hetero, Mg_hetero_se, 'Color', CBP(3, :), 'LineWidth', 2);
plot(Mg, Mg_hetero, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
xlabel('Cation concentration (mM)'); 
ylabel('Average FRET efficiency, \langleE\rangle (%)');
xticks([1e-3 10e-3 100e-3 1]); xticklabels({'1', '10', '100', '1000'});
axis([0.5e-3 4 -5 75]);

subplot(122)
set(gca, 'FontSize', 12, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, control_homo, control_homo_se, 'Color', CBP(5, :), ...
    'LineWidth', 2);
plot(Mg, control_homo, 'o', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(5, :));
errorbar(Mg, long_homo, long_homo_se, 'Color', CBP(7, :), 'LineWidth', 2);
plot(Mg, long_homo, 'd', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
errorbar(Mg, untethered_homo, untethered_homo_se, 'Color', CBP(7, :), ...
    'LineWidth', 2);
plot(Mg, untethered_homo, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
errorbar(Na, Na_homo, Na_homo_se, 'Color', CBP(2, :), 'LineWidth', 2);
plot(Na, Na_homo, '<', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(2, :));
errorbar(Ca, Ca_homo, Ca_homo_se, 'Color', CBP(6, :), 'LineWidth', 2);
plot(Ca, Ca_homo, 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(6, :));
errorbar(Mg, Mg_homo, Mg_homo_se, 'Color', CBP(7, :), 'LineWidth', 2);
plot(Mg, Mg_homo, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
xlabel('Cation concentration (mM)'); 
ylabel('Average FRET efficiency, \langleE\rangle (%)');
xticks([1e-3 10e-3 100e-3 1]); xticklabels({'1', '10', '100', '1000'});
axis([0.5e-3 4 -5 75]);