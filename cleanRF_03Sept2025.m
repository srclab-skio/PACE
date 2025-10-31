%% Random Forest with Second Derivatives, Validated on In-Situ Data
clear; close all; clc;

% Hardcode wavelengths from column names (nm)
wavelengths = [403,405,408,410,413,415,418,420,423,425,428,430,433,435,438,440,443,445,448,450,...
               453,455,458,460,463,465,468,470,473,475,478,480,483,485,488,490,493,495,498,500,...
               503,505,508,510,513,515,518,520,523,525,528,530,533,535,538,540,543,545,548,550,...
               553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
               625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
               675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];

num_wl = length(wavelengths);

% Load simulated data (training)
data_table_sim = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve'); 
% data_table_sim = readtable('Filtered_Subset_global_GLORIA_LAB_close_tss20.csv', 'VariableNamingRule', 'preserve'); 

tss_sim = data_table_sim.TSS;  % TSS values (mg/L)
rrs_data_sim = table2array(data_table_sim(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
num_samples_sim = size(rrs_data_sim, 1);

% Convert response variable to double and check for validity
tss_sim = double(tss_sim);
% Check for NaN or Inf in response variable
nan_idx_sim = isnan(tss_sim) | isinf(tss_sim);
if any(nan_idx_sim)
    fprintf('Found %d samples with NaN or Inf in TSS (simulated). Removing them.\n', sum(nan_idx_sim));
    tss_sim = tss_sim(~nan_idx_sim);
    rrs_data_sim = rrs_data_sim(~nan_idx_sim, :);
    num_samples_sim = size(rrs_data_sim, 1);
end
% Ensure Rrs is positive
rrs_data_sim = max(rrs_data_sim, eps);

% Load in-situ data (testing)
% data_table_insitu = readtable('Subset_Simulated_close_tss_rrs_spectra.csv', 'VariableNamingRule', 'preserve');
data_table_insitu = readtable('Filtered_Subset_global_GLORIA_LAB_close_tss5.csv', 'VariableNamingRule', 'preserve');
% data_table_insitu = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
% data_table_insitu = readtable('Cruise_and_Sapelo_trip_4_SRF_rrs.csv', 'VariableNamingRule', 'preserve');
tss_insitu = data_table_insitu.TSS;  % TSS values (mg/L)
rrs_data_insitu = table2array(data_table_insitu(:, 23:end));  % Rrs spectra
num_samples_insitu = size(rrs_data_insitu, 1);

% Convert response variable to double and check for validity
tss_insitu = double(tss_insitu);
% Check for NaN or Inf in response variable
nan_idx_insitu = isnan(tss_insitu) | isinf(tss_insitu);
if any(nan_idx_insitu)
    fprintf('Found %d samples with NaN or Inf in TSS (in-situ). Removing them.\n', sum(nan_idx_insitu));
    tss_insitu = tss_insitu(~nan_idx_insitu);
    rrs_data_insitu = rrs_data_insitu(~nan_idx_insitu, :);
    num_samples_insitu = size(rrs_data_insitu, 1);
end
% Ensure Rrs is positive
rrs_data_insitu = max(rrs_data_insitu, eps);

% Interpolate in-situ data to model wavelengths (assuming in-situ wavelengths are 400-800 nm)
insitu_wl = 400:800;  % Placeholder; replace with actual wavelengths if available
if size(rrs_data_insitu, 2) ~= num_wl
    fprintf('In-situ data has %d wavelengths; interpolating to match model wavelengths (403-718 nm).\n', size(rrs_data_insitu, 2));
    rrs_data_insitu_interp = zeros(num_samples_insitu, num_wl);
    for i = 1:num_samples_insitu
        rrs_data_insitu_interp(i, :) = interp1(insitu_wl, rrs_data_insitu(i, :), wavelengths, 'linear', 'extrap');
    end
    rrs_data_insitu = rrs_data_insitu_interp;
end

% Compute second derivatives for simulated data (training)
second_deriv_sim = zeros(size(rrs_data_sim));
for i = 1:num_samples_sim
    spec = rrs_data_sim(i, :);
    spec = sgolayfilt(spec, 3, 11); 
    % spec = smoothdata(spec, 'movmean', 5);  % Fallback
%     try
%         spec = sgolayfilt(spec, 3, 5);  % 3rd-order, 11-point window
%     catch
%         spec = smoothdata(spec, 'movmean', 5);  % Fallback
%     end
    first_deriv = gradient(spec, wavelengths);
    second_deriv_sim(i, :) = gradient(first_deriv, wavelengths);
end

% Check for NaN in second derivative features (simulated)
nan_in_features_sim = any(any(isnan(second_deriv_sim)));
fprintf('NaN in second derivative features (simulated): %d\n', nan_in_features_sim);
nan_samples_sim = any(isnan(second_deriv_sim), 2);
if any(nan_samples_sim)
    fprintf('Found %d samples with NaN in second derivatives (simulated). Removing them.\n', sum(nan_samples_sim));
    second_deriv_sim = second_deriv_sim(~nan_samples_sim, :);
    tss_sim = tss_sim(~nan_samples_sim);
    num_samples_sim = size(second_deriv_sim, 1);
end

% Compute second derivatives for in-situ data (testing)
second_deriv_insitu = zeros(size(rrs_data_insitu));
for i = 1:num_samples_insitu
    spec = rrs_data_insitu(i, :);
    try
        spec = sgolayfilt(spec, 3, 11);
    catch
        spec = smoothdata(spec, 'movmean', 5);
    end
    first_deriv = gradient(spec, wavelengths);
    second_deriv_insitu(i, :) = gradient(first_deriv, wavelengths);
end

% Check for NaN in second derivative features (in-situ)
nan_in_features_insitu = any(any(isnan(second_deriv_insitu)));
fprintf('NaN in second derivative features (in-situ): %d\n', nan_in_features_insitu);
nan_samples_insitu = any(isnan(second_deriv_insitu), 2);
if any(nan_samples_insitu)
    fprintf('Found %d samples with NaN in second derivatives (in-situ). Removing them.\n', sum(nan_samples_insitu));
    second_deriv_insitu = second_deriv_insitu(~nan_samples_insitu, :);
    tss_insitu = tss_insitu(~nan_samples_insitu);
    num_samples_insitu = size(second_deriv_insitu, 1);
end

% Use best hyperparameters from original code
best_rs = 49;
best_leaf = 6;
best_trees = 100;

% Train model on full simulated data
rng(best_rs);
rf_model = TreeBagger(best_trees, second_deriv_sim, tss_sim, ...
    'Method', 'regression', 'Options', statset('UseParallel', false),...
    'MinLeafSize', best_leaf,...
    'OOBPredictorImportance', 'on');
% 'Options', statset('UseParallel', false), for stopping parallel
% processing. But. it might take much time.
% Save the model
% save('PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_Smoothing.mat', ...
%      'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% save('PACE_rf_2nd_der_insitu20_Sim_validated.mat', ...
%      'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');

% Predict on in-situ test data
y_pred = predict(rf_model, second_deriv_insitu);

% Evaluate performance
eval_metrics = @(y_true, y_pred) struct(...
    'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
    'RMSE', sqrt(mean((y_true - y_pred).^2)));
metrics = eval_metrics(tss_insitu, y_pred);

% Additional metrics for high TSS (>50 mg/L)
high_tss_idx = tss_insitu > 50;
metrics_high = eval_metrics(tss_insitu(high_tss_idx), y_pred(high_tss_idx));

% Display metrics
fprintf('Random Forest Model (2nd Derivatives, In-Situ Validation):\n');
fprintf('Overall R²: %.4f, RMSE: %.4f mg/L\n', metrics.R2, metrics.RMSE);
fprintf('High TSS (>50 mg/L) R²: %.4f, RMSE: %.4f mg/L, n: %d\n', ...
        metrics_high.R2, metrics_high.RMSE, sum(high_tss_idx));

% Check test set diagnostics
fprintf('Test set size (in-situ): %d\n', length(tss_insitu));
fprintf('Any NaN in y_test: %d\n', any(isnan(tss_insitu)));
fprintf('Any NaN in y_pred: %d\n', any(isnan(y_pred)));

% Feature importances for best model
importances = rf_model.OOBPermutedPredictorDeltaError;
[sorted_imp, imp_idx] = sort(importances, 'descend');

fprintf('\nTop 50 Wavelengths by Importance (Random Forest, 2nd Derivatives):\n');
for i = 1:min(50, length(wavelengths))
    fprintf('Wavelength %d nm: Importance %.4f\n', wavelengths(imp_idx(i)), sorted_imp(i));
end

% Create feature importance plot
figure('Position', [100, 100, 1200, 600]);

% Subplot 1: Bar plot of top 50 important features
subplot(1, 2, 1);
barh(sorted_imp(1:min(50, length(wavelengths))), 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'YTick', 1:min(50, length(wavelengths)));
set(gca, 'YTickLabel', wavelengths(imp_idx(1:min(50, length(wavelengths)))));
ylabel('Wavelength (nm)');
xlabel('Feature Importance');
title('Top 50 Important Wavelengths (2nd Derivatives)');
grid on;
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');

% Subplot 2: Spectral plot of all feature importances
subplot(1, 2, 2);
plot(wavelengths, importances, 'b-', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
hold on;
[~, top20_idx] = maxk(importances, 20);
scatter(wavelengths(top20_idx), importances(top20_idx), 100, 'ro', 'filled', ...
        'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Feature Importance');
title('Feature Importance Spectrum (2nd Derivatives)');
grid on;
legend('All features', 'Top 20 features', 'Location', 'best');
for i = 1:20
    text(wavelengths(top20_idx(i)), importances(top20_idx(i)) + 0.1, ...
         sprintf('%d nm', wavelengths(top20_idx(i))), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
ylim_vals = [min(importances)*0.9, max(importances)*1.1];
subplot(1, 2, 1); xlim(ylim_vals);
subplot(1, 2, 2); ylim(ylim_vals);
sgtitle('Random Forest Feature Importance - 2nd Derivatives (In-Situ Validation)', ...
        'FontSize', 14, 'FontWeight', 'bold');
% saveas(gcf, 'feature_imp_2nd_der_insitu20_sim_validated.png');
filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_insitu20_Sim_valid" + best_rs + "_" + best_trees  + ".png");
exportgraphics(gcf,filename,'Resolution',600)
% filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_LAB5_" + best_rs + "_" + best_trees  + ".eps");
% exportgraphics(gcf,filename2);


% % Display feature importance summary
% fprintf('\nFeature Importance Summary:\n');
% fprintf('Mean importance: %.4f\n', mean(importances));
% fprintf('Max importance: %.4f (at %d nm)\n', max(importances), wavelengths(imp_idx(1)));
% fprintf('Min importance: %.4f (at %d nm)\n', min(importances), wavelengths(imp_idx(end)));
% fprintf('Importance range: %.4f\n', range(importances));

% % Model diagnostics
% fprintf('\nModel Diagnostics:\n');
% fprintf('Number of Observations (simulated, training): %d\n', num_samples_sim);
% fprintf('Number of Observations (in-situ, testing): %d\n', num_samples_insitu);

%% Simple Scatter plot for predictions
figure(15);
hold on;
scatter(tss_insitu, y_pred, 20, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
plot([min(tss_insitu), max(tss_insitu)], [min(tss_insitu), max(tss_insitu)], 'b--', 'LineWidth', 2);
r2 = 1 - sum((tss_insitu - y_pred).^2) / sum((tss_insitu - mean(tss_insitu)).^2);
rmse = sqrt(mean((tss_insitu - y_pred).^2));
combined_str = {sprintf('R² = %.3f', r2), ...
                sprintf('RMSE = %.2f mg/L', rmse), ...
                sprintf('n = %d', length(tss_insitu)), ...
                sprintf('Bias = %.2f mg/L', mean(y_pred - tss_insitu))};
text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
xlabel('In-Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-2, 100]);
ylim([-2, 100]);
title('In-situ vs. RF Predicted TSS', ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);
p = polyfit(tss_insitu, y_pred, 1);
regress_line = polyval(p, xlim);
plot(xlim, regress_line, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');
slope = p(1);
intercept = p(2);
legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope, intercept)}, ...
       'Location', 'southeast');
hold off;
set(gcf, 'Position', [100, 100, 800, 700]);
%saveas(gcf, 'scatter_rf_2nd_deriv_validated.png');
filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_LAB5_" + best_rs + "_" + best_trees  + ".png");
filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_LAB5_" + best_rs + "_" + best_trees  + ".eps");
%exportgraphics(gcf,filename,'Resolution',600)
%exportgraphics(gcf,filename2);

%%
%% Scatter plot for predictions with distinct colors and histograms
% Create scatter plot figure
figure(15);
hold on;

% Get unique countries and assign colors from parula colormap
countries = unique(data_table_insitu.Country);
num_countries = length(countries);
colors = parula(num_countries); % Use parula for distinct colors

% Plot scatter points for each country with 50% transparency
for i = 1:num_countries
    country_mask = strcmp(data_table_insitu.Country, countries{i});
    scatter(tss_insitu(country_mask), y_pred(country_mask), 20, colors(i,:), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
        'DisplayName', countries{i}, 'AlphaData', 0.5); % 50% transparency
end

% Define axis limits
x_limits = [-2, 100];
y_limits = [-2, 100];

% Plot 1:1 line across the entire figure
plot(x_limits, y_limits, 'b--', 'LineWidth', 2, 'DisplayName', '1:1 Line');

% Calculate regression line
p = polyfit(tss_insitu, y_pred, 1);
slope = p(1);
intercept = p(2);
regress_line = polyval(p, x_limits);
plot(x_limits, regress_line, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Regression: y = %.2fx + %.2f', slope, intercept));

% Calculate metrics
r2 = 1 - sum((tss_insitu - y_pred).^2) / sum((tss_insitu - mean(tss_insitu)).^2);
rmse = sqrt(mean((tss_insitu - y_pred).^2));
bias = mean(y_pred - tss_insitu);

% Create text box with metrics
combined_str = {sprintf('R² = %.3f', r2), ...
                sprintf('RMSE = %.2f mg/L', rmse), ...
                sprintf('n = %d', length(tss_insitu)), ...
                sprintf('Bias = %.2f mg/L', bias)};
text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

% Set labels and title for scatter plot
xlabel('In-Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
title('In-situ vs. RF Predicted TSS by Country', ...
      'FontSize', 14, 'FontWeight', 'bold');

% Set axis limits and grid
xlim(x_limits);
ylim(y_limits);
grid on;
set(gca, 'FontSize', 12);

% Add legend with 2 columns
legend('show', 'Location', 'southeast', 'NumColumns', 2);

% Add histogram of residuals on the left middle side
% residuals = tss_insitu - y_pred;
axes('Position', [0.15, 0.35, 0.15, 0.25]); % [left, bottom, width, height]
histogram(tss_insitu, 20, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
xlabel('Residuals (mg L^{-1})', 'FontSize', 10);
ylabel('Frequency', 'FontSize', 10);
set(gca, 'FontSize', 10);
grid on;

% Finalize scatter plot figure
hold off;
set(gcf, 'Position', [100, 100, 800, 700]);

% Save scatter plot figure
filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
    "testing_GLO_LAB5_" + best_rs + "_" + best_trees + ".png");
filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
    "testing_GLO_LAB5_" + best_rs + "_" + best_trees + ".eps");
exportgraphics(gcf, filename, 'Resolution', 600);
exportgraphics(gcf, filename2);

% Create separate figure for histograms
figure(16);
set(gcf, 'Position', [900, 100, 800, 400]);

% Subplot 1: Histogram of insitu
subplot(1, 3, 1);
histogram(tss_insitu, 20, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
xlabel('In situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 12, 'FontWeight', 'bold');
title('Histogram of in situ TSS', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);

% Subplot 2: Histogram of predicted TSS
subplot(1, 3, 2);
histogram(y_pred, 20, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
xlabel('Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 12, 'FontWeight', 'bold');
title('Histogram of Predicted TSS', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);

% Subplot 2: Histogram of Simulated TSS
subplot(1, 3, 3);
histogram(tss_sim, 20, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
xlabel('Simulated TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 12, 'FontWeight', 'bold');
title('Histogram of Simulated TSS', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);

% Save histogram figure
hist_filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
    "histograms_GLO_LAB5_" + best_rs + "_" + best_trees + ".png");
hist_filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
    "histograms_GLO_LAB5_" + best_rs + "_" + best_trees + ".eps");
exportgraphics(gcf, hist_filename, 'Resolution', 600);
exportgraphics(gcf, hist_filename2);

% %% Scatter plot for predictions with different colors for countries
% figure();
% hold on;
% 
% % Get unique countries and assign colors
% countries = unique(data_table_insitu.Country);
% num_countries = length(countries);
% colors = lines(num_countries); % Generate distinct colors for each country
% 
% % Plot scatter points for each country
% for i = 1:num_countries
%     country_mask = strcmp(data_table_insitu.Country, countries{i});
%     scatter(tss_insitu(country_mask), y_pred(country_mask), 20, colors(i,:), ...
%         'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
%         'DisplayName', countries{i});
% end
% 
% % Define axis limits
% x_limits = [-2, 100];
% y_limits = [-2, 100];
% 
% % Plot 1:1 line across the entire figure
% plot(x_limits, y_limits, 'b--', 'LineWidth', 2, 'DisplayName', '1:1 Line');
% 
% % Calculate regression line
% p = polyfit(tss_insitu, y_pred, 1);
% slope = p(1);
% intercept = p(2);
% regress_line = polyval(p, x_limits);
% plot(x_limits, regress_line, 'k:', 'LineWidth', 1.5, ...
%     'DisplayName', sprintf('Regression: y = %.2fx + %.2f', slope, intercept));
% 
% % Calculate metrics
% r2 = 1 - sum((tss_insitu - y_pred).^2) / sum((tss_insitu - mean(tss_insitu)).^2);
% rmse = sqrt(mean((tss_insitu - y_pred).^2));
% bias = mean(y_pred - tss_insitu);
% 
% % Create text box with metrics
% combined_str = {sprintf('R² = %.3f', r2), ...
%                 sprintf('RMSE = %.2f mg/L', rmse), ...
%                 sprintf('n = %d', length(tss_insitu)), ...
%                 sprintf('Bias = %.2f mg/L', bias)};
% text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
%     'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
%     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
% 
% % Set labels and title
% xlabel('In-Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
% title('In-situ vs. RF Predicted TSS by Country', ...
%       'FontSize', 14, 'FontWeight', 'bold');
% 
% % Set axis limits and grid
% xlim(x_limits);
% ylim(y_limits);
% grid on;
% set(gca, 'FontSize', 12);
% 
% % Add legend
% legend('show', 'Location', 'southeast', 'NumColumns', 2);
% 
% % Finalize figure
% hold off;
% set(gcf, 'Position', [100, 100, 800, 700]);
% 
% % Save figure
% filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
%     "Countries_testing_GLO_LAB5_" + best_rs + "_" + best_trees + ".png");
% filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
%     "Countries_testing_GLO_LAB5_" + best_rs + "_" + best_trees + ".eps");
% exportgraphics(gcf, filename, 'Resolution', 600);
% exportgraphics(gcf, filename2);
%%
%% Scatter plot for predictions with different colors for countries and histogram
% figure(15);
% hold on;
% 
% % Get unique countries and assign colors
% countries = unique(data_table_insitu.Country);
% num_countries = length(countries);
% colors = lines(num_countries); % Generate distinct colors for each country
% 
% % Plot scatter points for each country
% for i = 1:num_countries
%     country_mask = strcmp(data_table_insitu.Country, countries{i});
%     scatter(tss_insitu(country_mask), y_pred(country_mask), 20, colors(i,:), ...
%         'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
%         'DisplayName', countries{i});
% end
% 
% % Define axis limits
% x_limits = [-2, 100];
% y_limits = [-2, 100];
% 
% % Plot 1:1 line across the entire figure
% plot(x_limits, y_limits, 'b--', 'LineWidth', 2, 'DisplayName', '1:1 Line');
% 
% % Calculate regression line
% p = polyfit(tss_insitu, y_pred, 1);
% slope = p(1);
% intercept = p(2);
% regress_line = polyval(p, x_limits);
% plot(x_limits, regress_line, 'k:', 'LineWidth', 1.5, ...
%     'DisplayName', sprintf('Regression: y = %.2fx + %.2f', slope, intercept));
% 
% % Calculate metrics
% r2 = 1 - sum((tss_insitu - y_pred).^2) / sum((tss_insitu - mean(tss_insitu)).^2);
% rmse = sqrt(mean((tss_insitu - y_pred).^2));
% bias = mean(y_pred - tss_insitu);
% 
% % Create text box with metrics
% combined_str = {sprintf('R² = %.3f', r2), ...
%                 sprintf('RMSE = %.2f mg/L', rmse), ...
%                 sprintf('n = %d', length(tss_insitu)), ...
%                 sprintf('Bias = %.2f mg/L', bias)};
% text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
%     'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
%     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
% 
% % Set labels and title
% xlabel('In-Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
% title('In-situ vs. RF Predicted TSS by Country', ...
%       'FontSize', 14, 'FontWeight', 'bold');
% 
% % Set axis limits and grid
% xlim(x_limits);
% ylim(y_limits);
% grid on;
% set(gca, 'FontSize', 12);
% 
% % Add legend with 2 columns
% legend('show', 'Location', 'southeast', 'NumColumns', 2);
% 
% % Add histogram of residuals on the left middle side
% residuals = tss_insitu - y_pred;
% axes('Position', [0.15, 0.35, 0.15, 0.25]); % [left, bottom, width, height]
% histogram(tss_insitu, 20, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
% xlabel('Residuals (mg L^{-1})', 'FontSize', 10);
% ylabel('Frequency', 'FontSize', 10);
% set(gca, 'FontSize', 10);
% grid on;
% 
% % Finalize figure
% hold off;
% set(gcf, 'Position', [100, 100, 800, 700]);
% 
% % Save figure
% filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
%     "testing_GLO_LAB5_" + best_rs + "_" + best_trees + ".png");
% filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
%     "testing_GLO_LAB5_" + best_rs + "_" + best_trees + ".eps");
% exportgraphics(gcf, filename, 'Resolution', 600);
% exportgraphics(gcf, filename2);

%% Scatter plot for predictions; advanced, publication
% figure(15);
% hold on;
% 
% % Check if we have enough data points for density calculation
% if length(tss_insitu) < 2
%     % Not enough points for density calculation, use simple coloring
%     scatter_handle = scatter(tss_insitu, y_pred, 40, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
%     cbar_label = '';
%     use_density = false;
% else
%     % Use kde2d for proper density calculation
%     try
%         [density, X, Y] = kde2d([tss_insitu, y_pred], min(50, length(tss_insitu)));
%         
%         % Check if we have valid density data
%         if all(isnan(density(:))) || numel(unique(tss_insitu)) < 2 || numel(unique(y_pred)) < 2
%             error('Insufficient data for density calculation');
%         end
%         
%         % Interpolate density to each data point
%         density_at_points = interp2(X, Y, density, tss_insitu, y_pred, 'linear', 0);
%         
%         % Scatter plot with proper density-based coloring
%         scatter_handle = scatter(tss_insitu, y_pred, 40, density_at_points, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
%         colormap(parula); % Use parula colormap (no red)
%         
%         cbar_label = 'Point Density';
%         use_density = true;
%         
%     catch
%         % Fallback: simple distance-based coloring
%         dist_from_center = sqrt((tss_insitu - mean(tss_insitu)).^2 + (y_pred - mean(y_pred)).^2);
%         density_at_points = 1./(1 + dist_from_center); % Inverse distance as density proxy
%         
%         scatter_handle = scatter(tss_insitu, y_pred, 40, density_at_points, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
%         colormap(parula); % Use parula colormap (no red)
%         
%         cbar_label = 'Data Intensity';
%         use_density = true;
%     end
% end
% 
% % Set axis limits to ensure 1:1 line goes corner to corner
% axis_limits = [-2, 100, -2, 100];
% xlim(axis_limits(1:2));
% ylim(axis_limits(3:4));
% 
% % 1:1 line - RED DASHED from corner to corner
% plot(axis_limits(1:2), axis_limits(1:2), 'r--', 'LineWidth', 2.5);
% 
% % Calculate regression line (only if we have enough points)
% if length(tss_insitu) >= 2
%     [p, S] = polyfit(tss_insitu, y_pred, 1);
%     regress_x = linspace(axis_limits(1), axis_limits(2), 100);
%     regress_y = polyval(p, regress_x);
%     
%     % Regression line - BLACK DASHED line
%     plot(regress_x, regress_y, 'k--', 'LineWidth', 2.5);
%     
%     % 95% confidence intervals - GREEN lines
%     residuals = y_pred - polyval(p, tss_insitu);
%     ci_range = 1.96 * std(residuals);
%     upper_ci = regress_y + ci_range;
%     lower_ci = regress_y - ci_range;
%     plot(regress_x, upper_ci, ':', 'LineWidth', 1.8, 'Color', [0, 0.5, 0]);
%     plot(regress_x, lower_ci, ':', 'LineWidth', 1.8, 'Color', [0, 0.5, 0]);
% else
%     p = [0, 0]; % Default values if not enough data
% end
% 
% % Calculate metrics
% if ~isempty(tss_insitu) && length(tss_insitu) >= 2
%     r2 = 1 - sum((tss_insitu - y_pred).^2) / sum((tss_insitu - mean(tss_insitu)).^2);
%     rmse = sqrt(mean((tss_insitu - y_pred).^2));
%     bias = mean(y_pred - tss_insitu);
%     n_points = length(tss_insitu);
% else
%     r2 = 0;
%     rmse = 0;
%     bias = 0;
%     n_points = length(tss_insitu);
% end
% 
% % SINGLE TEXT BOX with all metrics - positioned between x=10 to 40, y=95 to 80
% combined_str = {sprintf('R² = %.2f', r2), ...
%                 sprintf('RMSE = %.2f mg/L', rmse), ...
%                 sprintf('n = %d', n_points), ...
%                 sprintf('Bias = %.2f mg/L', bias)};
%             
% text_box = annotation('textbox', [0.2, 0.5, 0.7, 0.7], 'String', combined_str, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
%     'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
%     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FitBoxToText','on');
% 
% % Axis labels with bold and thick formatting
% xlabel('In-Situ TSS (mg L^{-1})', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 14, 'FontWeight', 'bold');
% 
% % Title
% title('In-situ vs. RF Predicted TSS', 'FontSize', 16, 'FontWeight', 'bold');
% 
% % Remove grid
% grid off;
% 
% % Make a square box with thick lines on all four sides
% set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2.5, ...
%          'Box', 'on', 'XColor', [0, 0, 0], 'YColor', [0, 0, 0]);
% 
% % Add colorbar only if we have density data
% if use_density && ~isempty(cbar_label)
%     cbar = colorbar('Location', 'east');
%     cbar.Label.String = cbar_label;
%     cbar.Label.FontSize = 12;
%     cbar.Label.FontWeight = 'bold';
%     cbar.FontWeight = 'bold';
%     
%     % Position colorbar between x=80 to 95 (right side of plot)
%     current_axis_pos = get(gca, 'Position');
%     cbar_width = 0.02;
%     cbar_height = 0.3;
%     cbar_x = current_axis_pos(1) + current_axis_pos(3) + 0.015; % Right side + small offset
%     cbar_y = current_axis_pos(2) + (current_axis_pos(4) - cbar_height)/2; % Vertically centered
%     set(cbar, 'Position', [cbar_x, cbar_y, cbar_width, cbar_height]);
% end
% 
% % Add regression equation at TOP-MIDDLE of the plot (only if we have regression)
% if length(tss_insitu) >= 2
%     text(50, 95, sprintf('y = %.2fx + %.2f', p(1), p(2)), ...
%          'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white', ...
%          'HorizontalAlignment', 'center', 'EdgeColor', 'black');
% end
% 
% % Legend (only add entries that exist)
% legend_entries = {};
% if ~isempty(tss_insitu)
%     legend_entries{end+1} = 'Data Points';
% end
% legend_entries{end+1} = '1:1 Line';
% if length(tss_insitu) >= 2
%     legend_entries{end+1} = 'Regression Line';
%     legend_entries{end+1} = '95% CI';
% end
% 
% if ~isempty(legend_entries)
%     legend(legend_entries, 'Location', 'southeast', 'FontSize', 10, 'Box', 'on');
% end
% 
% hold off;
% 
% % Ensure equal aspect ratio for proper square shape
% axis equal;
% axis([-2, 100, -2, 100]);
% 
% % Figure position and saving
% set(gcf, 'Position', [100, 100, 800, 700]);
% % saveas(gcf, 'scatter_rf_2nd_deriv_sapelo_cruise_matchups.png');
% %%
% filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_LAB5_" + best_rs + "_" + best_trees  + ".png");
% filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_LAB5_" + best_rs + "_" + best_trees  + ".eps");
% exportgraphics(gcf,filename,'Resolution',600)
% exportgraphics(gcf,filename2);

%% Different geographical area: 31 OCT 2025
%% Random Forest with Second Derivatives, Validated on In-Situ Data
clear; close all; clc;

% Hardcode wavelengths from column names (nm)
wavelengths = [403,405,408,410,413,415,418,420,423,425,428,430,433,435,438,440,443,445,448,450,...
               453,455,458,460,463,465,468,470,473,475,478,480,483,485,488,490,493,495,498,500,...
               503,505,508,510,513,515,518,520,523,525,528,530,533,535,538,540,543,545,548,550,...
               553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
               625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
               675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];

num_wl = length(wavelengths);

% Load simulated data (training)
data_table_sim = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve'); 
tss_sim = data_table_sim.TSS;  % TSS values (mg/L)
rrs_data_sim = table2array(data_table_sim(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
num_samples_sim = size(rrs_data_sim, 1);

% Convert response variable to double and check for validity
tss_sim = double(tss_sim);
nan_idx_sim = isnan(tss_sim) | isinf(tss_sim);
if any(nan_idx_sim)
    fprintf('Found %d samples with NaN or Inf in TSS (simulated). Removing them.\n', sum(nan_idx_sim));
    tss_sim = tss_sim(~nan_idx_sim);
    rrs_data_sim = rrs_data_sim(~nan_idx_sim, :);
    num_samples_sim = size(rrs_data_sim, 1);
end
% Ensure Rrs is positive
rrs_data_sim = max(rrs_data_sim, eps);

% Load in-situ data (testing)
data_table_insitu = readtable('Filtered_Subset_global_GLORIA_LAB_close_tss5.csv', 'VariableNamingRule', 'preserve');
tss_insitu = data_table_insitu.TSS;  % TSS values (mg/L)
rrs_data_insitu = table2array(data_table_insitu(:, 23:end));  % Rrs spectra
latitude_full = data_table_insitu.Latitude;
longitude_full = data_table_insitu.Longitude;
country_full = data_table_insitu.Country;  
num_samples_insitu = size(rrs_data_insitu, 1);

% Convert response variable to double and check for validity
tss_insitu = double(tss_insitu);
nan_idx_insitu = isnan(tss_insitu) | isinf(tss_insitu);
if any(nan_idx_insitu)
    fprintf('Found %d samples with NaN or Inf in TSS (in-situ). Removing them.\n', sum(nan_idx_insitu));
    tss_insitu = tss_insitu(~nan_idx_insitu);
    rrs_data_insitu = rrs_data_insitu(~nan_idx_insitu, :);
    latitude_full = latitude_full(~nan_idx_insitu);
    longitude_full = longitude_full(~nan_idx_insitu);
    country_full = country_full(~nan_idx_insitu);
    num_samples_insitu = size(rrs_data_insitu, 1);
end
% Ensure Rrs is positive
rrs_data_insitu = max(rrs_data_insitu, eps);

% Interpolate in-situ data to model wavelengths (assuming in-situ wavelengths are 400-800 nm)
insitu_wl = 400:800;
if size(rrs_data_insitu, 2) ~= num_wl
    fprintf('In-situ data has %d wavelengths; interpolating to match model wavelengths (403-718 nm).\n', size(rrs_data_insitu, 2));
    rrs_data_insitu_interp = zeros(num_samples_insitu, num_wl);
    for i = 1:num_samples_insitu
        rrs_data_insitu_interp(i, :) = interp1(insitu_wl, rrs_data_insitu(i, :), wavelengths, 'linear', 'extrap');
    end
    rrs_data_insitu = rrs_data_insitu_interp;
end

% Compute second derivatives for simulated data (training)
second_deriv_sim = zeros(size(rrs_data_sim));
for i = 1:num_samples_sim
    spec = rrs_data_sim(i, :);
    spec = sgolayfilt(spec, 3, 11);
    first_deriv = gradient(spec, wavelengths);
    second_deriv_sim(i, :) = gradient(first_deriv, wavelengths);
end

% Check for NaN in second derivative features (simulated)
nan_in_features_sim = any(any(isnan(second_deriv_sim)));
fprintf('NaN in second derivative features (simulated): %d\n', nan_in_features_sim);
nan_samples_sim = any(isnan(second_deriv_sim), 2);
if any(nan_samples_sim)
    fprintf('Found %d samples with NaN in second derivatives (simulated). Removing them.\n', sum(nan_samples_sim));
    second_deriv_sim = second_deriv_sim(~nan_samples_sim, :);
    tss_sim = tss_sim(~nan_samples_sim);
    num_samples_sim = size(second_deriv_sim, 1);
end

% Compute second derivatives for in-situ data (testing)
second_deriv_insitu = zeros(size(rrs_data_insitu));
for i = 1:num_samples_insitu
    spec = rrs_data_insitu(i, :);
    try
        spec = sgolayfilt(spec, 3, 11);
    catch
        spec = smoothdata(spec, 'movmean', 5);
    end
    first_deriv = gradient(spec, wavelengths);
    second_deriv_insitu(i, :) = gradient(first_deriv, wavelengths);
end

% Check for NaN in second derivative features (in-situ)
nan_in_features_insitu = any(any(isnan(second_deriv_insitu)));
fprintf('NaN in second derivative features (in-situ): %d\n', nan_in_features_insitu);
nan_samples_insitu = any(isnan(second_deriv_insitu), 2);
if any(nan_samples_insitu)
    fprintf('Found %d samples with NaN in second derivatives (in-situ). Removing them.\n', sum(nan_samples_insitu));
    second_deriv_insitu = second_deriv_insitu(~nan_samples_insitu, :);
    tss_insitu = tss_insitu(~nan_samples_insitu);
    latitude_full = latitude_full(~nan_samples_insitu);
    longitude_full = longitude_full(~nan_samples_insitu);
    country_full = country_full(~nan_samples_insitu);
    num_samples_insitu = size(second_deriv_insitu, 1);
end

% Use best hyperparameters from original code
best_rs = 49;
best_leaf = 6;
best_trees = 100;

% Train model on full simulated data
rng(best_rs);
rf_model = TreeBagger(best_trees, second_deriv_sim, tss_sim, ...
    'Method', 'regression', 'Options', statset('UseParallel', false),...
    'MinLeafSize', best_leaf,...
    'OOBPredictorImportance', 'on');
% Save the model
save('PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat', ...
     'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% Predict on in-situ test data
y_pred = predict(rf_model, second_deriv_insitu);

% Evaluate performance
eval_metrics = @(y_true, y_pred) struct(...
    'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
    'RMSE', sqrt(mean((y_true - y_pred).^2)));
metrics = eval_metrics(tss_insitu, y_pred);

% Additional metrics for high TSS (>50 mg/L)
high_tss_idx = tss_insitu > 50;
metrics_high = eval_metrics(tss_insitu(high_tss_idx), y_pred(high_tss_idx));

% Display metrics
fprintf('Random Forest Model (2nd Derivatives, In-Situ Validation):\n');
fprintf('Overall R²: %.4f, RMSE: %.4f mg/L\n', metrics.R2, metrics.RMSE);
fprintf('High TSS (>50 mg/L) R²: %.4f, RMSE: %.4f mg/L, n: %d\n', ...
        metrics_high.R2, metrics_high.RMSE, sum(high_tss_idx));

% Check test set diagnostics
fprintf('Test set size (in-situ): %d\n', length(tss_insitu));
fprintf('Any NaN in y_test: %d\n', any(isnan(tss_insitu)));
fprintf('Any NaN in y_pred: %d\n', any(isnan(y_pred)));

%% Feature importances for best model
importances = rf_model.OOBPermutedPredictorDeltaError;
[sorted_imp, imp_idx] = sort(importances, 'descend');

fprintf('\nTop 50 Wavelengths by Importance (Random Forest, 2nd Derivatives):\n');
for i = 1:min(50, length(wavelengths))
    fprintf('Wavelength %d nm: Importance %.4f\n', wavelengths(imp_idx(i)), sorted_imp(i));
end

% Create feature importance plot
figure('Position', [100, 100, 1200, 600]);
subplot(1, 2, 1);
barh(sorted_imp(1:min(50, length(wavelengths))), 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'YTick', 1:min(50, length(wavelengths)));
set(gca, 'YTickLabel', wavelengths(imp_idx(1:min(50, length(wavelengths)))));
ylabel('Wavelength (nm)');
xlabel('Feature Importance');
title('Top 50 Important Wavelengths (2nd Derivatives)');
grid on;
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');

subplot(1, 2, 2);
plot(wavelengths, importances, 'b-', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
hold on;
[~, top20_idx] = maxk(importances, 20);
scatter(wavelengths(top20_idx), importances(top20_idx), 100, 'ro', 'filled', ...
        'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Feature Importance');
title('Feature Importance Spectrum (2nd Derivatives)');
grid on;
legend('All features', 'Top 20 features', 'Location', 'best');
for i = 1:20
    text(wavelengths(top20_idx(i)), importances(top20_idx(i)) + 0.1, ...
         sprintf('%d nm', wavelengths(top20_idx(i))), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
ylim_vals = [min(importances)*0.9, max(importances)*1.1];
subplot(1, 2, 1); xlim(ylim_vals);
subplot(1, 2, 2); ylim(ylim_vals);
sgtitle('Random Forest Feature Importance - 2nd Derivatives (In-Situ Validation)', ...
        'FontSize', 14, 'FontWeight', 'bold');
filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_insitu20_Sim_valid" + best_rs + "_" + best_trees  + ".png");
drawnow;
exportgraphics(gcf, filename, 'Resolution', 600);
close(gcf);

%% Simple Scatter plot for predictions (overall)
figure(15);
hold on;
scatter(tss_insitu, y_pred, 20, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
plot([min(tss_insitu), max(tss_insitu)], [min(tss_insitu), max(tss_insitu)], 'b--', 'LineWidth', 2);
r2 = 1 - sum((tss_insitu - y_pred).^2) / sum((tss_insitu - mean(tss_insitu)).^2);
rmse = sqrt(mean((tss_insitu - y_pred).^2));
combined_str = {sprintf('R² = %.3f', r2), ...
                sprintf('RMSE = %.2f mg/L', rmse), ...
                sprintf('n = %d', length(tss_insitu)), ...
                sprintf('Bias = %.2f mg/L', mean(y_pred - tss_insitu))};
text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
xlabel('In-Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-2, 100]);
ylim([-2, 100]);
title('In-situ vs. RF Predicted TSS', ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);
p = polyfit(tss_insitu, y_pred, 1);
regress_line = polyval(p, xlim);
plot(xlim, regress_line, 'k:', 'LineWidth', 1.5);
slope = p(1);
intercept = p(2);
legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope, intercept)}, ...
       'Location', 'southeast');
hold off;
set(gcf, 'Position', [100, 100, 800, 700]);
%%
filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_LAB5_" + best_rs + "_" + best_trees  + ".png");
filename2 = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', "testing_GLO_LAB5_" + best_rs + "_" + best_trees  + ".eps");
drawnow;
%exportgraphics(gcf,filename,'Resolution',600)
%exportgraphics(gcf,filename2);
close(gcf);

%% Categorize in-situ data based on latitude and longitude
regions = cell(num_samples_insitu, 1);
for i = 1:num_samples_insitu
    lat = latitude_full(i);
    lon = longitude_full(i);
    
    if lat >= 25 && lat <= 35 && lon >= -84 && lon <= -78
        regions{i} = 'South Atlantic Bight';
    elseif lon >= -85 && lon <= -65 && lat >= 25 && lat <= 45
        regions{i} = 'US East Coast';
    elseif strcmp(country_full{i}, 'USA')
        regions{i} = 'Rest of the US';
    elseif abs(lat) <= 10
        regions{i} = 'Equatorial Regions';
    elseif lat > 0
        if lon >= -130 && lon <= -115 && lat >= 30 && lat <= 50
            regions{i} = 'US West Coast';
        elseif lon >= 120 && lon <= 150 && lat >= 30 && lat <= 45
            regions{i} = 'North Pacific';
        elseif lon >= -10 && lon <= 30 && lat >= 35 && lat <= 60
            regions{i} = 'Europe';
        else
            regions{i} = 'Other Northern Hemisphere';
        end
    else  % lat < 0
        if lon >= 110 && lon <= 160 && lat >= -40 && lat <= -10
            regions{i} = 'Australia';
        elseif lon >= 165 && lon <= 180 && lat >= -50 && lat <= -30
            regions{i} = 'New Zealand';
        elseif lon >= -60 && lon <= -30 && lat >= -30 && lat <= -5
            regions{i} = 'Brazil';
        elseif lon >= 15 && lon <= 35 && lat >= -35 && lat <= -20
            regions{i} = 'South Africa';
        elseif lon <= -100 || lon >= 100
            regions{i} = 'South Pacific';
        else
            regions{i} = 'Other Southern Hemisphere';
        end
    end
end

% Scatter plots for different geographical regions as subplots
unique_regions = unique(regions);
num_regions = length(unique_regions);
rows = 3;
cols = 4;

figure(16);
set(gcf, 'Position', [100, 100, 1200, 800]);

% Determine global axis limits
global_min = min([min(tss_insitu), min(y_pred)]);
global_max = max([max(tss_insitu), max(y_pred)]);
global_limits = [global_min - 2, global_max + 2];

% Create tight subplot layout using axes
margin = 0.02; % Small margin for axes
width = (1 - margin * (cols + 1)) / cols;
height = (1 - margin * (rows + 1)) / rows;

for r = 1:min(num_regions, rows * cols)
    reg = unique_regions{r};
    idx = strcmp(regions, reg);
    if sum(idx) < 1
        continue;  % Skip if no samples in this region
    end
    
    tss_reg = tss_insitu(idx);
    pred_reg = y_pred(idx);
    
    % Calculate axes position [left bottom width height]
    col_idx = mod(r-1, cols);
    row_idx = floor((r-1)/cols);
    left = margin + col_idx * (width + margin);
    bottom = 1 - (row_idx + 1) * (height + margin);
    ax = axes('Position', [left, bottom, width, height]);
    
    hold on;
    scatter(tss_reg, pred_reg, 20, 'b', 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    
    plot(global_limits, global_limits, 'b--', 'LineWidth', 1.5);
    
    r2 = 1 - sum((tss_reg - pred_reg).^2) / sum((tss_reg - mean(tss_reg)).^2);
    rmse = sqrt(mean((tss_reg - pred_reg).^2));
    bias = mean(pred_reg - tss_reg);
    combined_str = sprintf('%s\nR² = %.2f\nRMSE = %.2f mg/L\nn = %d\nBias = %.2f mg/L', reg, r2, rmse, sum(idx), bias);
    text(0.05, 0.95, combined_str, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
         'FontSize', 8, 'FontWeight', 'bold');
    
    xlabel('In-Situ TSS (mg L^{-1})', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 10, 'FontWeight', 'bold');
    xlim(global_limits);
    ylim(global_limits);
    set(gca, 'FontSize', 8, 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');
    axis square;
    
    p = polyfit(tss_reg, pred_reg, 1);
    regress_line = polyval(p, global_limits);
    plot(global_limits, regress_line, 'k:', 'LineWidth', 1.5);
    slope = p(1);
    intercept = p(2);
    legend({'Data Points', '1:1 Line', sprintf('y = %.2fx + %.2f', slope, intercept)}, ...
           'Location', 'southeast', 'FontSize', 7);
    hold off;
end

%sgtitle('In-situ vs. RF Predicted TSS by Region', 'FontSize', 14, 'FontWeight', 'bold');
%% Broad category 31 OCT 2025
%% Categorize in-situ data into EXACTLY 4 redefined regions
regions = cell(num_samples_insitu, 1);

% ----- Core region definitions (for fallback) -----
% [min_lon  max_lon  min_lat  max_lat  region_name]
core_regions = {
    -130, -65,  25,  70, 'North America';       % USA/Canada + coasts
    -85,   45, -60, -20, 'Southern Hemisphere'; % South America, S. Africa
    110,  160, -50, -10, 'Southern Hemisphere'; % Australia
    165,  180, -50, -30, 'Southern Hemisphere'; % New Zealand
    -10,   40,  35,  70, 'European Regions';    % Europe + UK + Norway
     90,  150,   0,  45, 'North Pacific Regions'
};

for i = 1:num_samples_insitu
    lat = latitude_full(i);
    lon = longitude_full(i);
    country = country_full{i};
    
    % ----- Normalise longitude to [-180 180] -----
    if lon > 180,  lon = lon - 360;  end
    
    assigned = false;
    
    % ----- 1. North America (country or coastal boxes) -----
    if strcmp(country,'USA') || strcmp(country,'CAN')
        regions{i} = 'North America';
        assigned = true;
    elseif (lat>=25 && lat<=35 && lon>=-84 && lon<=-78) || ... % South Atlantic Bight
           (lat>=25 && lat<=45 && lon>=-85 && lon<=-65) || ... % US East Coast
           (lat>=30 && lat<=50 && lon>=-130 && lon<=-115)      % US West Coast
        regions{i} = 'North America';
        assigned = true;
    end
    
    % ----- 2. European Regions -----
    if ~assigned && lon>=-10 && lon<=40 && lat>=35 && lat<=70
        regions{i} = 'European Regions';
        assigned = true;
    end
    
    % ----- 3. North Pacific Regions -----
    if ~assigned && lon>=90 && lon<=150 && lat>=0 && lat<=45
        regions{i} = 'North Pacific Regions';
        assigned = true;
    end
    
    % ----- 4. Southern Hemisphere landmasses -----
    if ~assigned && lat<0
        if (lon>=110 && lon<=160 && lat>=-50 && lat<=-10) || ... % Australia
           (lon>=165 && lon<=180 && lat>=-50 && lat<=-30) || ... % New Zealand
           (lon>=-85 && lon<=-30 && lat>=-60 && lat<=-5)  || ... % South America
           (lon>=10  && lon<=45  && lat>=-35 && lat<=-20)        % South Africa
            regions{i} = 'Southern Hemisphere';
            assigned = true;
        end
    end
    
    % ----- Fallback: nearest core region centre -----
    if ~assigned
        min_dist = inf;  best_region = '';
        for k = 1:size(core_regions,1)
            clon = (core_regions{k,1}+core_regions{k,2})/2;
            clat = (core_regions{k,3}+core_regions{k,4})/2;
            dist = sqrt((lon-clon)^2 + (lat-clat)^2);
            if dist < min_dist
                min_dist = dist;  best_region = core_regions{k,5};
            end
        end
        regions{i} = best_region;
    end
end

%% ---------- Scatter plots (4 regions, 2x2 layout) ----------
unique_regions = unique(regions);               % should be exactly 4
rows = 2;  cols = 2;
figure(16);  clf;
set(gcf,'Position',[100 100 1000 800]);

% Global axis limits
global_min = min([tss_insitu(:); y_pred(:)]);
global_max = max([tss_insitu(:); y_pred(:)]);
global_limits = [max(0,global_min-2) , global_max+2];

% ----- Region order & colours (RGB triplets) -----
plot_order = {'North America','European Regions','North Pacific Regions','Southern Hemisphere'};
hex_colors = {'#1f77b4','#ff7f0e','#2ca02c','#d62728'};   % same as before
rgb_colors = cellfun(@(c)sscanf(c(2:end),'%2x%2x%2x')'/255, hex_colors, 'UniformOutput',false);
rgb_colors = vertcat(rgb_colors{:});                     % 4x3 matrix

margin = 0.06;
width  = (1 - margin*(cols+1))/cols;
height = (1 - margin*(rows+1))/rows;

for r = 1:numel(plot_order)
    reg = plot_order{r};
    idx = strcmp(regions,reg);
    if ~any(idx), continue; end
    
    tss_reg = tss_insitu(idx);
    pred_reg = y_pred(idx);
    
    % ----- Axes position -----
    col_idx = mod(r-1,cols);
    row_idx = floor((r-1)/cols);
    left   = margin + col_idx*(width + margin);
    bottom = 1 - (row_idx+1)*(height + margin);
    ax = axes('Position',[left bottom width height]);  hold on;
    
    % ----- Scatter (RGB colour) -----
    scatter(tss_reg, pred_reg, 30, rgb_colors(r,:), 'filled', ...
            'MarkerFaceAlpha',0.7,'MarkerEdgeColor','k','LineWidth',0.5);
    
    % 1:1 line
    plot(global_limits,global_limits,'k--','LineWidth',1.5);
    
    % Regression line
    p = polyfit(tss_reg,pred_reg,1);
    plot(global_limits,polyval(p,global_limits),'r-','LineWidth',1.8);
    
    % Statistics
    r2   = 1 - sum((tss_reg-pred_reg).^2)/sum((tss_reg-mean(tss_reg)).^2);
    rmse = sqrt(mean((tss_reg-pred_reg).^2));
    bias = mean(pred_reg-tss_reg);
    n    = numel(tss_reg);
    
    % Annotation
    txt = sprintf('%s\nR^2 = %.3f\nRMSE = %.2f\nn = %d\nBias = %.2f',...
                  reg,r2,rmse,n,bias);
    text(0.05,0.95,txt,'Units','normalized','VerticalAlignment','top',...
         'HorizontalAlignment','left','BackgroundColor','w','EdgeColor','k',...
         'FontSize',8,'FontWeight','bold','Margin',3);
    
    % Labels / formatting
    xlabel('In-Situ TSS (mg L^{-1})','FontSize',10);
    ylabel('Predicted TSS (mg L^{-1})','FontSize',10);
    title(reg,'FontSize',11,'FontWeight','bold','Color',rgb_colors(r,:));
    xlim(global_limits); ylim(global_limits);
    set(gca,'FontSize',9,'Box','on','LineWidth',0.8);
    axis square; grid off;
    
    % Legend
    legend({'Data','1:1 Line',sprintf('Fit (y = %.2fx + %.2f)',p(1),p(2))},...
           'Location','southeast','FontSize',7,'Box','off');
    hold off;
end

sgtitle('In-situ vs. RF Predicted TSS by 4 Global Regions',...
        'FontSize',14,'FontWeight','bold');
tightfig(gcf);
%%  Save figure
filename = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
                    sprintf("4_scatter_regions_%d_%d.png", best_rs, best_trees));
drawnow;
exportgraphics(gcf, filename, 'Resolution', 600);
filename_eps = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
                        sprintf("4_scatter_regions_%d_%d.eps", best_rs, best_trees));
drawnow;
exportgraphics(gcf, filename_eps);
close(gcf);

% %% Global map plot
% figure(17);
% set(gcf, 'Position', [100, 100, 1200, 800]);
% ga = geoaxes('Basemap', 'grayland');
% hold on;
% 
% % Color by country using a discrete colormap
% unique_countries = unique(country_full);
% num_countries = length(unique_countries);
% colors = hsv(num_countries); % Use lines, hsv, parula, turbo for discrete color scale
% 
% for c = 1:num_countries
%     cnt = unique_countries{c};
%     idx = strcmp(country_full, cnt);
%     geoscatter(latitude_full(idx), longitude_full(idx), 20, colors(c,:), 'filled');
% end
% 
% title('Global Distribution of In-Situ Data Points by Country', 'FontSize', 14, 'FontWeight', 'bold');
% legend(unique_countries, 'Location', 'northeast', 'NumColumns', 3, 'FontSize', 10);
% hold off;
% 
% %%  Save figure
% filename_map = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
%                         sprintf("global_regional_map_%d_%d.png", best_rs, best_trees));
% drawnow;
% exportgraphics(gcf, filename_map, 'Resolution', 600);
% filename_map_eps = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
%                             sprintf("global_regional_map_%d_%d.eps", best_rs, best_trees));
% drawnow;
% exportgraphics(gcf, filename_map_eps);
% close(gcf);

%%
%% --------------------------------------------------------------
%% Global map – colour by the 4 regions (same colours as scatter plots)
%% --------------------------------------------------------------
figure(17);  clf;
set(gcf,'Position',[100 100 1300 800],'Color','w');

%--- Geoaxes ---------------------------------------------------------
ga = geoaxes('Basemap','grayland');
hold(ga,'on');

%--- SAME 4 colours (RGB triplets) ------------------------------------
hex = {'#1f77b4','#ff7f0e','#2ca02c','#d62728'};      % identical to scatter plots
rgb = cellfun(@(c)sscanf(c(2:end),'%2x%2x%2x')'/255, hex, ...
              'UniformOutput',false);
rgb = vertcat(rgb{:});                               % 4×3 matrix

%--- Fixed region order (must match scatter-plot order) ---------------
plot_order = {'North America','European Regions',...
              'North Pacific Regions','Southern Hemisphere'};

%--- Plot points ------------------------------------------------------
legend_handles = gobjects(0);
for r = 1:numel(plot_order)
    reg = plot_order{r};
    idx = strcmp(regions,reg);
    if ~any(idx), continue; end
    
    h = geoscatter(ga, latitude_full(idx), longitude_full(idx), ...
                   40, rgb(r,:), 'filled', ...
                   'MarkerEdgeColor','k','LineWidth',1.2, ...
                   'MarkerFaceAlpha',0.8);
    legend_handles(end+1) = h;   %#ok<AGROW>
end

%--- Title (bold & large) --------------------------------------------
title(ga,'Global Distribution of In-Situ Data Points by Region',...
      'FontSize',16,'FontWeight','bold','Color','k');

% %--- Grid (R2022a-compatible) ----------------------------------------
% ga.Grid = 'on';                     % turn grid on
% ga.GridLineWidth = 1.2;             % thick grid lines
% ga.GridColor = [0.3 0.3 0.3];        % dark grey
% ga.MinorGrid = 'off';

%--- Axis labels (bold) ----------------------------------------------
ga.LongitudeLabel.FontWeight = 'bold';
ga.LatitudeLabel.FontWeight  = 'bold';
ga.LongitudeLabel.FontSize   = 11;
ga.LatitudeLabel.FontSize    = 11;

%--- Full view --------------------------------------------------------
geolimits(ga,[-70 70],[-180 180]);

%--- Legend (bold, thick box, outside) -------------------------------
legend(ga,legend_handles,plot_order,...
       'Location','north','Orientation','horizontal',...
       'NumColumns',4,'FontSize',11,'FontWeight','bold',...
       'Box','on','EdgeColor','k','LineWidth',1.5);

%--- Final map styling ------------------------------------------------
geobasemap(ga,'grayland');
set(ga,'LineWidth',1.5);   % thick coastlines
hold(ga,'off');

% %% Recompute regions if not already available (from your first code snippet)
% %Uncomment and adjust if 'regions' is not precomputed
% num_samples_insitu = length(latitude_full);  % Adjust if needed
% regions = cell(num_samples_insitu, 1);
% for i = 1:num_samples_insitu
%     lat = latitude_full(i);
%     lon = longitude_full(i);
%     
%     if lat >= 25 && lat <= 35 && lon >= -84 && lon <= -78
%         regions{i} = 'South Atlantic Bight';
%     elseif lon >= -85 && lon <= -65 && lat >= 25 && lat <= 45
%         regions{i} = 'US East Coast';
%     elseif strcmp(country_full{i}, 'USA')  % Assumes country_full is available
%         regions{i} = 'Rest of the US';
%     elseif abs(lat) <= 10
%         regions{i} = 'Equatorial Regions';
%     elseif lat > 0
%         if lon >= -130 && lon <= -115 && lat >= 30 && lat <= 50
%             regions{i} = 'US West Coast';
%         elseif lon >= 120 && lon <= 150 && lat >= 30 && lat <= 45
%             regions{i} = 'North Pacific';
%         elseif lon >= -10 && lon <= 30 && lat >= 35 && lat <= 60
%             regions{i} = 'Europe';
%         else
%             regions{i} = 'Other Northern Hemisphere';
%         end
%     else  % lat < 0
%         if lon >= 110 && lon <= 160 && lat >= -40 && lat <= -10
%             regions{i} = 'Australia';
%         elseif lon >= 165 && lon <= 180 && lat >= -50 && lat <= -30
%             regions{i} = 'New Zealand';
%         elseif lon >= -60 && lon <= -30 && lat >= -30 && lat <= -5
%             regions{i} = 'Brazil';
%         elseif lon >= 15 && lon <= 35 && lat >= -35 && lat <= -20
%             regions{i} = 'South Africa';
%         elseif lon <= -100 || lon >= 100
%             regions{i} = 'South Pacific';
%         else
%             regions{i} = 'Other Southern Hemisphere';
%         end
%     end
% end
% 
% %% Global map: Color by geographical region (using the same categorization)
% figure(17);
% set(gcf, 'Position', [100, 100, 1200, 800]);
% ga = geoaxes('Basemap', 'grayland');
% hold on;
% 
% % Color by region using a discrete colormap
% unique_regions = unique(regions);
% num_regions = length(unique_regions);
% colors = hsv(num_regions);  % Discrete colormap (swap to lines, hsv, parula, turbo if preferred)
% 
% for r = 1:num_regions
%     reg = unique_regions{r};
%     idx = strcmp(regions, reg);
%     if sum(idx) < 1
%         continue;  % Skip if no samples in this region
%     end
%     geoscatter(latitude_full(idx), longitude_full(idx), 20, colors(r,:), 'filled');
% end
% 
% title('Global Distribution of In-Situ Data Points by Region', 'FontSize', 14, 'FontWeight', 'bold');
% legend(unique_regions, 'Location', 'northeast', 'NumColumns', 3, 'FontSize', 10);
% hold off;

%% Save figure
filename_map = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
                        sprintf("4global_regional_map_by_region_%d_%d.png", best_rs, best_trees));
drawnow;
exportgraphics(gcf, filename_map, 'Resolution', 600);
filename_map_eps = fullfile('/Users/masud/OneDriveUGA/QWIP/Figs', ...
                            sprintf("4global_regional_map_by_region_%d_%d.eps", best_rs, best_trees));
drawnow;
exportgraphics(gcf, filename_map_eps);
% close(gcf);

%% Find samples close to 1:1 line (differences <10)

% close_idx = abs(y_pred - tss_insitu) < 10;
% num_close = sum(close_idx);
% fprintf('Found %d samples with TSS difference <10 mg/L.\n', num_close);
% 
% % Write close samples to CSV
% if num_close > 0
%     close_table = data_table_insitu(close_idx, :);
%     close_table.Predicted_TSS = y_pred(close_idx);
%     writetable(close_table, 'Subset_global_GLORIA_LAB_close_tss.csv');
%     fprintf('Saved close samples with all variables and predicted TSS to Subset_global_GLORIA_LAB_close_tss.csv\n');
% else
%     fprintf('No samples with TSS difference <10 mg/L found.\n');
% end

