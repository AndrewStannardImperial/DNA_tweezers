% Andrew Stannard
clear variables; close all;

CBP = [0 114 187; 241 90 35; 0 168 117; 235 222 58; ...
    219 111 171; 248 147 29; 1 186 242; 189 190 192]/255;

Mg = [1 2 5 10 20 50 100 200];
Mg_F_mean_DA_hetero = [218058 214644 210761 202770 191970 167039 136392 ...
    112864];
Mg_F_se_DA_hetero = [3495 3587 3961 3858 4304 6328 7461 8052];
Mg_F_mean_DO_hetero = [215885 215079 210902 204138 201203 194804 185947 ...
    179816];
Mg_F_se_DO_hetero = [4109 3625 5587 3439 2014 1869 3840 2703];
Mg_F_mean_DA_homo = [218875 212285 209379 200676 184169 146206 111582 ...
    86028];
Mg_F_se_DA_homo = [368 1942 1615 2012 2883 2720 2891 2972];
Mg_F_mean_DO_homo = [218282 214914 211939 206603 201017 191619 184452 ...
    175872];
Mg_F_se_DO_homo = [2440 3010 3185 2935 3034 3036 2498 1423];

Ca = [1 2 5 10 20 50 100 200];
Ca_F_mean_DA_hetero = [218572 211197 203984 196579 173230 136341 99327 ...
    79091];
Ca_F_se_DA_hetero = [2194 2581 4279 2982 5031 3202 3140 7359];
Ca_F_mean_DO_hetero = [216339 212688 205362 200037 190685 184408 171954 ...
    165977];
Ca_F_se_DO_hetero = [1316 770 1638 909 1044 385 3847 1999];
Ca_F_mean_DA_homo = [215040 209160 199521 186803 159749 108785 75401 ...
    56696];
Ca_F_se_DA_homo = [4059 3829 3652 2798 2492 3520 3703 3720];
Ca_F_mean_DO_homo = [214777 210928 203627 197211 189504 181470 172373 ...
    164899];
Ca_F_se_DO_homo = [3473 3855 3909 3814 3360 3185 3388 2453];

subplot(231)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, Mg_F_mean_DA_hetero, Mg_F_se_DA_hetero, ...
    'Color', CBP(3, :), 'LineWidth', 2);
plot(Mg, Mg_F_mean_DA_hetero, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
errorbar(Mg, Mg_F_mean_DO_hetero, Mg_F_se_DO_hetero, ...
    'Color', CBP(3, :), 'LineWidth', 2);
plot(Mg, Mg_F_mean_DO_hetero, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(3, :));
errorbar(Mg, Mg_F_mean_DA_homo, Mg_F_se_DA_homo, ...
    'Color', CBP(7, :), 'LineWidth', 2);
plot(Mg, Mg_F_mean_DA_homo, '^', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
errorbar(Mg, Mg_F_mean_DO_homo, Mg_F_se_DO_homo, ...
    'Color', CBP(7, :), 'LineWidth', 2);
plot(Mg, Mg_F_mean_DO_homo, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(7, :));
ylabel('Fluorescence (a.u.)'); xlabel('[Mg^{2+}] (mM)');
xticks([1 10 100]); xticklabels({'1', '10', '100'});
xlim([0.5 400]); ylim([36197 237545]);

subplot(232)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Ca, Ca_F_mean_DA_hetero, Ca_F_se_DA_hetero, ...
    'Color', CBP(4, :), 'LineWidth', 2);
plot(Ca, Ca_F_mean_DA_hetero, 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(4, :));
errorbar(Ca, Ca_F_mean_DO_hetero, Ca_F_se_DO_hetero, ...
    'Color', CBP(4, :), 'LineWidth', 2);
plot(Ca, Ca_F_mean_DO_hetero, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(4, :));
errorbar(Ca, Ca_F_mean_DA_homo, Ca_F_se_DA_homo, ...
    'Color', CBP(6, :), 'LineWidth', 2);
plot(Ca, Ca_F_mean_DA_homo, 'v', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(6, :));
errorbar(Ca, Ca_F_mean_DO_homo, Ca_F_se_DO_homo, ...
    'Color', CBP(6, :), 'LineWidth', 2);
plot(Ca, Ca_F_mean_DO_homo, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', CBP(6, :));
ylabel('Fluorescence (a.u.)'); xlabel('[Ca^{2+}] (mM)');
xticks([1 10 100]); xticklabels({'1', '10', '100'});
xlim([0.5 400]); ylim([36197 237545]);

Mg_F_mean = [217084 214997 211421 205370 201110 193212 185199 177844];
Mg_F_se = [2000 1885 2585 1907 1458 1626 1868 1558];

Ca_F_mean = [215558 211808 204495 198624 190094 182939 172163 165438];
Ca_F_se = [1534 1631 1748 1714 1437 1471 2053 1293];

subplot(233)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(Mg, Mg_F_mean, Mg_F_se, ...
    'Color', (CBP(3, :) + CBP(7, :)) / 2, 'LineWidth', 2);
plot(Mg, Mg_F_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(3, :) + CBP(7, :)) / 2);
errorbar(Ca, Ca_F_mean, Ca_F_se, ...
    'Color', (CBP(4, :) + CBP(6, :)) / 2, 'LineWidth', 2);
plot(Ca, Ca_F_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(4, :) + CBP(6, :)) / 2);
ylabel('DO_{pooled} fluorescence (a.u.)'); xlabel('[M^{2+}] (mM)');
xticks([1 10 100]); xticklabels({'1', '10', '100'});
xlim([0.5 400]); ylim([36197 237545]);

Mg_corr_mean = 0.972; Mg_corr_se = 0.008;
Mg_RQY_mean = Mg_F_mean/Mg_F_mean(1)*Mg_corr_mean;
Mg_RQY_se = Mg_RQY_mean.*sqrt((Mg_F_se(1)/Mg_F_mean(1))^2+...
    (Mg_corr_se/Mg_corr_mean)^2+(Mg_F_se./Mg_F_mean).^2);

Ca_corr_mean = 0.957; Ca_corr_se = 0.023;
Ca_RQY_mean = Ca_F_mean/Ca_F_mean(1)*Ca_corr_mean;
Ca_RQY_se = Ca_RQY_mean.*sqrt((Ca_F_se(1)/Ca_F_mean(1))^2+...
    (Ca_corr_se/Ca_corr_mean)^2+(Ca_F_se./Ca_F_mean).^2);

subplot(234)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
plot(0.25, 1, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', [0.75 0.75 0.75]);
errorbar(Mg, Mg_RQY_mean, Mg_RQY_se, ...
    'Color', (CBP(3, :) + CBP(7, :)) / 2, 'LineWidth', 2, 'LineStyle', ...
    'none');
plot(Mg, Mg_RQY_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(3, :) + CBP(7, :)) / 2);
errorbar(Ca, Ca_RQY_mean, Ca_RQY_se, ...
    'Color', (CBP(4, :) + CBP(6, :)) / 2, 'LineWidth', 2, 'LineStyle', ...
    'none');
plot(Ca, Ca_RQY_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(4, :) + CBP(6, :)) / 2);
ylabel('Relative quantum yield'); xlabel('[M^{2+}] (mM)');
xticks([0.1 1 10 100]); xticklabels({'0', '1', '10', '100'});
xlim([0.125 400]); ylim([0.6867 1.0285]);

Mg_RFR_mean = (Mg_RQY_mean).^(1/6);
Mg_RFR_se = Mg_RFR_mean.*Mg_RQY_se/6./Mg_RQY_mean;
Mg_x = logspace(log10(0.5),log10(400),101);
Mg_model = 1 - 8.056e-3*log(1+0.438*Mg_x);

Ca_RFR_mean = (Ca_RQY_mean).^(1/6);
Ca_RFR_se = Ca_RFR_mean.*Ca_RQY_se/6./Ca_RQY_mean;
Ca_x = logspace(log10(0.5),log10(400),101);
Ca_model = 1 - 9.626e-3*log(1+0.846*Mg_x);

subplot(235)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
plot(0.25, 1, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', [0.75 0.75 0.75]);
plot(Mg_x, Mg_model, 'Color', (CBP(3, :) + CBP(7, :)) / 2, 'LineWidth', 2);
errorbar(Mg, Mg_RFR_mean, Mg_RFR_se, ...
    'Color', (CBP(3, :) + CBP(7, :)) / 2, 'LineWidth', 2, 'LineStyle', ...
    'none');
plot(Mg, Mg_RFR_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(3, :) + CBP(7, :)) / 2);
plot(Ca_x, Ca_model, 'Color', (CBP(4, :) + CBP(6, :)) / 2, 'LineWidth', 2);
errorbar(Ca, Ca_RFR_mean, Ca_RFR_se, ...
    'Color', (CBP(4, :) + CBP(6, :)) / 2, 'LineWidth', 2, 'LineStyle', ...
    'none');
plot(Ca, Ca_RFR_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(4, :) + CBP(6, :)) / 2);
ylabel('Relative Forster radius'); xlabel('[M^{2+}] (mM)');
xticks([0.1 1 10 100]); xticklabels({'0', '1', '10', '100'});
xlim([0.125 400]); ylim([0.9403 1.0054]);

Mg_FR_mean = 5.71 * Mg_RFR_mean;
Mg_FR_se = Mg_FR_mean .* sqrt((Mg_RFR_se./Mg_RFR_mean).^2+(0.02/5.71)^2);

Ca_FR_mean = 5.71 * Ca_RFR_mean;
Ca_FR_se = Ca_FR_mean .* sqrt((Ca_RFR_se./Ca_RFR_mean).^2+(0.02/5.71)^2);

subplot(236)
set(gca, 'FontSize', 14, 'XScale', 'log'); hold on; axis square;
errorbar(0.25, 5.71, 0.02, ...
    'Color', [0.75 0.75 0.75], 'LineWidth', 2);
plot(Mg_x, 5.71*Mg_model, 'Color', (CBP(3, :) + CBP(7, :)) / 2, ...
    'LineWidth', 2);
plot(0.25, 5.71, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', [0.75 0.75 0.75]);
errorbar(Mg, Mg_FR_mean, Mg_FR_se, ...
    'Color', (CBP(3, :) + CBP(7, :)) / 2, 'LineWidth', 2, ...
    'LineStyle', 'none');
plot(Mg, Mg_FR_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(3, :) + CBP(7, :)) / 2);
plot(Ca_x, 5.71*Ca_model, 'Color', (CBP(4, :) + CBP(6, :)) / 2, ...
    'LineWidth', 2);
errorbar(Ca, Ca_FR_mean, Ca_FR_se, ...
    'Color', (CBP(4, :) + CBP(6, :)) / 2, 'LineWidth', 2, 'LineStyle', ...
    'none');
plot(Ca, Ca_FR_mean, 's', 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 10, 'MarkerFaceColor', (CBP(4, :) + CBP(6, :)) / 2);
ylabel('Forster radius (nm)'); xlabel('[M^{2+}] (mM)');
xticks([0.1 1 10 100]); xticklabels({'0', '1', '10', '100'});
xlim([0.125 400]); ylim([5.3597 5.7637]);