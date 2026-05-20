% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

load Fig3BC_fit_data.txt

x = Fig3BC_fit_data(:, 1);
Mg_fit = Fig3BC_fit_data(:, 2:10);
Ca_fit = Fig3BC_fit_data(:, 11:end);

conc = [1 2 5 10 20 50 100 200];
Mg_hetero_E = [-1.03 0.18 0.02 0.67 4.59 14.23 26.64 37.25]/100;
Mg_hetero_E_e = [1.16 1.58 1.47 0.84 1.65 2.95 3.13 3.43]/100;
Ca_hetero_E = [-1.04 0.70 0.64 1.72 9.17 26.06 42.14 52.27]/100;
Ca_hetero_E_e = [1.49 1.57 2.86 1.89 2.25 1.76 3.15 4.99]/100;
Mg_homo_E = [-0.29 1.21 1.19 2.86 8.38 23.70 39.52 51.10]/100;
Mg_homo_E_e = [1.03 0.61 0.78 0.57 0.82 0.78 0.85 1.39]/100;
Ca_homo_E = [-0.12 0.84 2.01 5.27 15.68 39.99 56.18 65.56]/100;
Ca_homo_E_e = [0.35 0.54 1.00 0.45 1.26 2.93 3.03 2.76]/100;

Ca_hetero_K = Ca_hetero_E./(1-Ca_hetero_E);
Ca_hetero_K_e = (Ca_hetero_K./Ca_hetero_E).^2.*Ca_hetero_E_e;
Mg_hetero_K = Mg_hetero_E./(1-Mg_hetero_E);
Mg_hetero_K_e = (Mg_hetero_K./Mg_hetero_E).^2.*Mg_hetero_E_e;
Ca_homo_K = Ca_homo_E./(1-Ca_homo_E);
Ca_homo_K_e = (Ca_homo_K./Ca_homo_E).^2.*Ca_homo_E_e;
Mg_homo_K = Mg_homo_E./(1-Mg_homo_E);
Mg_homo_K_e = (Mg_homo_K./Mg_homo_E).^2.*Mg_homo_E_e;

Ca_hetero_DG = -log(Ca_hetero_K);
Ca_hetero_DG_e = Ca_hetero_K_e./Ca_hetero_K;
Ca_homo_DG = -log(Ca_homo_K);
Ca_homo_DG_e = Ca_homo_K_e./Ca_homo_K;
Mg_hetero_DG = -log(Mg_hetero_K);
Mg_hetero_DG_e = Mg_hetero_K_e./Mg_hetero_K;
Mg_homo_DG = -log(Mg_homo_K);
Mg_homo_DG_e = Mg_homo_K_e./Mg_homo_K;

subplot(131);
set(gca, 'FontSize', 12, 'XScale', 'log','Ydir','reverse'); 
hold on; axis square;
h1 = fill([x; flipud(x)], [Mg_fit(:,4); flipud(Mg_fit(:,6))], ...
    CBP(3,:), 'LineStyle', 'none'); set(h1, 'facealpha', .25);
plot(x, Mg_fit(:,5), '--', 'LineWidth', 2, 'Color', CBP(3,:));
h2 = fill([x; flipud(x)], [Mg_fit(:,1); flipud(Mg_fit(:,3))], ...
    CBP(7,:), 'LineStyle', 'none'); set(h2, 'facealpha', .25);
plot(x, Mg_fit(:,2), '--', 'LineWidth', 2, 'Color', CBP(7,:));
errorbar(conc(5:end), Mg_hetero_DG(5:end), Mg_hetero_DG_e(5:end), ...
    'Color', CBP(3, :), 'LineStyle', 'none', 'LineWidth', 2);
plot(conc(5:end), Mg_hetero_DG(5:end), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
errorbar(conc(5:end), Mg_homo_DG(5:end), Mg_homo_DG_e(5:end), ...
    'Color', CBP(7, :), 'LineStyle', 'none', 'LineWidth', 2);
plot(conc(5:end), Mg_homo_DG(5:end), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
xlabel('[Mg^{2+}] (mM)'); 
ylabel('Coalignment free energy, \DeltaG (k_BT)');
xticks([20 50 100 200]); xticklabels({'20', '50', '100', '200'});
yticks([-1 0 1 2 3 4]); yticklabels({'-1', '0', '1', '2', '3', '4'});
xlim([10 400]); ylim([-1 4]);

subplot(132);
set(gca, 'FontSize', 12, 'XScale', 'log','Ydir','reverse'); 
hold on; axis square;
h3 = fill([x; flipud(x)], [Ca_fit(:,4); flipud(Ca_fit(:,6))], ...
    CBP(4,:), 'LineStyle', 'none'); set(h3, 'facealpha', .25);
plot(x, Ca_fit(:,5), '--', 'LineWidth', 2, 'Color', CBP(4,:));
h4 = fill([x; flipud(x)], [Ca_fit(:,1); flipud(Ca_fit(:,3))], ...
    CBP(6,:), 'LineStyle', 'none'); set(h4, 'facealpha', .25);
plot(x, Ca_fit(:,2), '--', 'LineWidth', 2, 'Color', CBP(6,:));
errorbar(conc(5:end), Ca_hetero_DG(5:end), Ca_hetero_DG_e(5:end), ...
    'Color', CBP(4, :), 'LineStyle', 'none', 'LineWidth', 2);
plot(conc(5:end), Ca_hetero_DG(5:end), 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(4, :));
errorbar(conc(5:end), Ca_homo_DG(5:end), Ca_homo_DG_e(5:end), ...
    'Color', CBP(6, :), 'LineStyle', 'none', 'LineWidth', 2);
plot(conc(5:end), Ca_homo_DG(5:end), 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(6, :));
xlabel('[Ca^{2+}] (mM)'); 
ylabel('Coalignment free energy, \DeltaG (k_BT)');
xticks([20 50 100 200]); xticklabels({'20', '50', '100', '200'});
yticks([-1 0 1 2 3 4]); yticklabels({'-1', '0', '1', '2', '3', '4'});
xlim([10 400]); ylim([-1 4]);

Ca_recog = -log(Ca_homo_K./Ca_hetero_K);
Ca_recog_e = sqrt((Ca_homo_K_e./Ca_homo_K).^2 + ...
    (Ca_hetero_K_e./Ca_hetero_K).^2);
Mg_recog = -log(Mg_homo_K./Mg_hetero_K);
Mg_recog_e = sqrt((Mg_homo_K_e./Mg_homo_K).^2 + ...
    (Mg_hetero_K_e./Mg_hetero_K).^2);

subplot(133);
set(gca, 'FontSize', 12, 'XScale', 'log','Ydir','reverse'); 
hold on; axis square;
h5 = fill([x; flipud(x)], [Ca_fit(:,7); flipud(Ca_fit(:,9))], ...
    CBP(6,:), 'LineStyle', 'none'); set(h5, 'facealpha', .25);
plot(x, Ca_fit(:,8), '--', 'LineWidth', 2, 'Color', CBP(6,:));
h6 = fill([x; flipud(x)], [Mg_fit(:,7); flipud(Mg_fit(:,9))], ...
    CBP(7,:), 'LineStyle', 'none'); set(h6, 'facealpha', .25);
plot(x, Mg_fit(:,8), '--', 'LineWidth', 2, 'Color', CBP(7,:));
errorbar(conc(5:end), Ca_recog(5:end), Ca_recog_e(5:end), 'Color', ...
    CBP(6, :), 'LineStyle', 'none', 'LineWidth', 2);
plot(conc(5:end), Ca_recog(5:end), 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(6, :));
errorbar(conc(5:end), Mg_recog(5:end), Mg_recog_e(5:end), ...
    'Color', CBP(7, :), 'LineStyle', 'none', 'LineWidth', 2);
plot(conc(5:end), Mg_recog(5:end), '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
xlabel('[M^{2+}] (mM)'); 
ylabel('Recognition free energy, \Delta\DeltaG_{recog} (k_BT)');
xticks([20 50 100 200]); xticklabels({'20', '50', '100', '200'});
yticks([-1.5 -1 -0.5 0 0.5]); 
yticklabels({'-1.5', '-1.0', '-0.5', '0.0', '0.5'});
xlim([10 400]); ylim([-1.5 0.5]); 