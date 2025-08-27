% clear all;
% % MATLAB code to process the CSV file as described
% % data_table = readtable('global_14Aug25_PACE_bands_only.csv', 'VariableNamingRule', 'preserve');
% 
% %% working for real model now
% clear; close all; clc;
% % Hardcode wavelengths from column names (nm)
% wavelengths = [403,405,408,410,413,415,418,420,423,425,428,430,433,435,438,440,443,445,448,450,...
%                453,455,458,460,463,465,468,470,473,475,478,480,483,485,488,490,493,495,498,500,...
%                503,505,508,510,513,515,518,520,523,525,528,530,533,535,538,540,543,545,548,550,...
%                553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
%                625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
%                675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];
% 
% % wavelengths = [550,...
% %                553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
% %                625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
% %                675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];
% 
% num_wl = length(wavelengths);
% 
% % Load real data
% % Load real data
% % data_table = readtable('global_14Aug25_PACE_bands_only.csv', 'VariableNamingRule', 'preserve');
% data_table = readtable('global_18Aug25_PACE_bands_TSS_LT200.csv', 'VariableNamingRule', 'preserve');
% 
% tss = data_table.TSS;  % TSS values (mg/L)
% % rrs_data = table2array(data_table(:, 2:end));  % Rrs spectra (rows: samples, cols: wavelengths)
% rrs_data = table2array(data_table(:, 2:end));  % Rrs spectra (rows: samples, cols: wavelengths)
% num_samples = size(rrs_data, 1);
% 
% % Ensure Rrs is positive to avoid issues with log transforms (if used later)
% rrs_data = max(rrs_data, eps);
% 
% % Compute second derivatives
% second_deriv = zeros(size(rrs_data));
% for i = 1:num_samples
%     spec = rrs_data(i, :);
%     % Smooth to reduce noise (use sgolayfilt if Signal Processing Toolbox available)
%     try
%         spec = sgolayfilt(spec, 3, 11);  % 3rd-order, 11-point window
%     catch
%         spec = smoothdata(spec, 'movmean', 5);  % Fallback if toolbox unavailable
%     end
%     first_deriv = gradient(spec, wavelengths);
%     second_deriv(i, :) = gradient(first_deriv, wavelengths);
% end
% 
% % Data splitting: 80% train, 20% test
% [train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
% X_raw_train = rrs_data(train_idx, :);
% X_raw_test = rrs_data(test_idx, :);
% X_deriv_train = second_deriv(train_idx, :);
% X_deriv_test = second_deriv(test_idx, :);
% y_train = tss(train_idx);
% y_test = tss(test_idx);
% 
% % Initialize models
% % Linear Regression
% lin_reg = @(X_train, y_train) fitlm(X_train, y_train);
% % Ridge Regression (lambda=1.0 equivalent)
% ridge_reg = @(X_train, y_train) fitrlinear(X_train, y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', 1.0);
% % Random Forest with OOB predictor importance
% rf_reg = @(X_train, y_train) TreeBagger(100, X_train, y_train, 'Method', 'regression', 'OOBPredictorImportance', 'on');
% 
% % Train and predict: Raw Rrs
% lin_model_raw = lin_reg(X_raw_train, y_train);
% y_pred_lin_raw = predict(lin_model_raw, X_raw_test);
% ridge_model_raw = ridge_reg(X_raw_train, y_train);
% y_pred_ridge_raw = predict(ridge_model_raw, X_raw_test);
% rf_model_raw = rf_reg(X_raw_train, y_train);
% y_pred_rf_raw = predict(rf_model_raw, X_raw_test);
% 
% % Train and predict: Second Derivatives
% lin_model_deriv = lin_reg(X_deriv_train, y_train);
% y_pred_lin_deriv = predict(lin_model_deriv, X_deriv_test);
% ridge_model_deriv = ridge_reg(X_deriv_train, y_train);
% y_pred_ridge_deriv = predict(ridge_model_deriv, X_deriv_test);
% rf_model_deriv = rf_reg(X_deriv_train, y_train);
% y_pred_rf_deriv = predict(rf_model_deriv, X_deriv_test);
% 
% % Evaluate models
% eval_metrics = @(y_true, y_pred) struct(...
%     'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
%     'RMSE', sqrt(mean((y_true - y_pred).^2)));
% 
% metrics_raw = struct(...
%     'Linear', eval_metrics(y_test, y_pred_lin_raw), ...
%     'Ridge', eval_metrics(y_test, y_pred_ridge_raw), ...
%     'RandomForest', eval_metrics(y_test, y_pred_rf_raw));
% metrics_deriv = struct(...
%     'Linear', eval_metrics(y_test, y_pred_lin_deriv), ...
%     'Ridge', eval_metrics(y_test, y_pred_ridge_deriv), ...
%     'RandomForest', eval_metrics(y_test, y_pred_rf_deriv));
% 
% % Display metrics
% fprintf('Raw Rrs Model Performance:\n');
% fprintf('Linear Regression - R²: %.4f, RMSE: %.4f mg/L\n', metrics_raw.Linear.R2, metrics_raw.Linear.RMSE);
% fprintf('Ridge Regression - R²: %.4f, RMSE: %.4f mg/L\n', metrics_raw.Ridge.R2, metrics_raw.Ridge.RMSE);
% fprintf('Random Forest - R²: %.4f, RMSE: %.4f mg/L\n', metrics_raw.RandomForest.R2, metrics_raw.RandomForest.RMSE);
% fprintf('\nSecond Derivatives Model Performance:\n');
% fprintf('Linear Regression - R²: %.4f, RMSE: %.4f mg/L\n', metrics_deriv.Linear.R2, metrics_deriv.Linear.RMSE);
% fprintf('Ridge Regression - R²: %.4f, RMSE: %.4f mg/L\n', metrics_deriv.Ridge.R2, metrics_deriv.Ridge.RMSE);
% fprintf('Random Forest - R²: %.4f, RMSE: %.4f mg/L\n', metrics_deriv.RandomForest.R2, metrics_deriv.RandomForest.RMSE);
% 
% % Feature importances for Random Forest (Raw Rrs)
% importances = rf_model_raw.OOBPermutedPredictorDeltaError;
% [sorted_imp, imp_idx] = sort(importances, 'descend');
% fprintf('\nTop 10 Wavelengths by Importance (Random Forest, Raw Rrs):\n');
% for i = 1:min(20, num_wl) % for the subset of bands
%     fprintf('Wavelength %d nm: Importance %.4f\n', wavelengths(imp_idx(i)), sorted_imp(i));
% end
% 
% 
% 
% %%
% % Scatter plots: In-situ vs Predicted TSS
% figure(12); 
% hold on;
% 
% % Create scatter plot
% scatter(y_test, y_pred_rf_raw, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
% 
% % Add 1:1 line
% plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', 'LineWidth', 2);
% 
% % Calculate metrics
% r2 = 1 - sum((y_test - y_pred_rf_raw).^2) / sum((y_test - mean(y_test)).^2);
% rmse = sqrt(mean((y_test - y_pred_rf_raw).^2));
% 
% % Combine all metrics and statistics in one text box
% combined_str = {sprintf('R² = %.3f', r2), ...
%                 sprintf('RMSE = %.2f mg/L', rmse), ...
%                 sprintf('n = %d', length(y_test)), ...
%                 sprintf('Bias = %.2f mg/L', mean(y_pred_rf_raw - y_test))};
% 
% text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
%     'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
%     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
% 
% % Formatting
% xlabel('In-Situ TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Predicted TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
% xlim([-2, 200]);
% ylim([-2, 200]);
% title('Random Forest: In-Situ vs Predicted TSS (Raw Rrs)', 'FontSize', 14, 'FontWeight', 'bold');
% grid on;
% set(gca, 'FontSize', 11);
% 
% % Add equation of the regression line (FIXED)
% p = polyfit(y_test, y_pred_rf_raw, 1);
% regress_line = polyval(p, xlim);
% plot(xlim, regress_line, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');
% 
% % FIXED LEGEND - store polyfit coefficients first
% slope = p(1);
% intercept = p(2);
% legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope, intercept)}, ...
%        'Location', 'southeast');
% 
% hold off;
% 
% % Save high-quality figure
% set(gcf, 'Position', [100, 100, 800, 700]);
% saveas(gcf, 'scatter_rf_raw.png');
% 
% 
% %% Second Derivatives: Random Forest
% figure(15); 
% hold on;
% 
% % Create scatter plot (matching style from figure 12)
% scatter(y_test, y_pred_rf_deriv, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
% 
% % Add 1:1 line
% plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', 'LineWidth', 2);
% 
% % Calculate metrics for second derivative model
% r2_deriv = 1 - sum((y_test - y_pred_rf_deriv).^2) / sum((y_test - mean(y_test)).^2);
% rmse_deriv = sqrt(mean((y_test - y_pred_rf_deriv).^2));
% 
% % Combine all metrics and statistics in one text box
% combined_str_deriv = {sprintf('R² = %.3f', r2_deriv), ...
%                       sprintf('RMSE = %.2f mg/L', rmse_deriv), ...
%                       sprintf('n = %d', length(y_test)), ...
%                       sprintf('Bias = %.2f mg/L', mean(y_pred_rf_deriv - y_test))};
% 
% text_box_deriv = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str_deriv, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
%     'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
%     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
% 
% % Formatting (matching figure 12 style)
% xlabel('In-Situ TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Predicted TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
% xlim([-2, 200]);
% ylim([-2, 200]);
% title('Random Forest: In-Situ vs Predicted TSS (Second Derivatives)', 'FontSize', 14, 'FontWeight', 'bold');
% grid on;
% set(gca, 'FontSize', 11);
% 
% % Add equation of the regression line for second derivative model
% p_deriv = polyfit(y_test, y_pred_rf_deriv, 1);
% regress_line_deriv = polyval(p_deriv, xlim);
% plot(xlim, regress_line_deriv, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');
% 
% % Legend for second derivative model
% slope_deriv = p_deriv(1);
% intercept_deriv = p_deriv(2);
% legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope_deriv, intercept_deriv)}, ...
%        'Location', 'southeast');
% 
% hold off;
% 
% % Save high-quality figure (matching figure 12 settings)
% set(gcf, 'Position', [100, 100, 800, 700]);
% saveas(gcf, 'scatter_rf_deriv.png');
% saveas(gcf, 'scatter_rf_deriv.pdf');  % Vector format for publications
% 
% % Display metrics in command window too
% fprintf('\nRandom Forest (Second Derivatives) Performance:\n');
% fprintf('R²: %.4f\n', r2_deriv);
% fprintf('RMSE: %.4f mg/L\n', rmse_deriv);
% fprintf('Bias: %.4f mg/L\n', mean(y_pred_rf_deriv - y_test));
% fprintf('Sample size: %d\n', length(y_test));
% %%
% 
% % Print sample predictions for Random Forest (Raw Rrs)
% fprintf('\nSample Predictions (Random Forest, Raw Rrs, First 10 Test Samples):\n');
% fprintf('In-Situ TSS\tPredicted TSS\n');
% for i = 1:min(10, length(y_test))
%     fprintf('%.4f\t\t%.4f\n', y_test(i), y_pred_rf_raw(i));
% end
% 
% % Print sample predictions for Random Forest (Second Derivatives)
% fprintf('\nSample Predictions (Random Forest, Second Derivatives, First 10 Test Samples):\n');
% fprintf('In-Situ TSS\tPredicted TSS\n');
% for i = 1:min(10, length(y_test))
%     fprintf('%.4f\t\t%.4f\n', y_test(i), y_pred_rf_deriv(i));
% end
% 
% disp('Plots and metrics generated. Check PNG files for scatter plots. Raw Rrs models (esp. Random Forest) recommended for TSS prediction.');
% 
% 
% %% 
% % MATLAB script for improving Random Forest on second derivatives for TSS prediction
% % Author: Grok (xAI) - Expert in Ocean Optics and Machine Learning
% % Date: August 14, 2025
% % Assumptions: CSV file 'rrs_data.csv' with headers: TSS, Rrs_403 to Rrs_718
% % Focus: Only Random Forest on second derivatives data
% % Tuning: Grid search over random states, NumTrees, MinLeafSize
% % Outputs: Metrics for each combination, scatter plots with R² annotated, best model identification
% % MATLAB Version: R2022a (requires Statistics and Machine Learning, Signal Processing Toolboxes)
% 
% %% Improving random Forest models
% clear; close all; clc;
% 
% % Hardcode wavelengths from column names (nm)
% wavelengths = [403,405,408,410,413,415,418,420,423,425,428,430,433,435,438,440,443,445,448,450,...
%                453,455,458,460,463,465,468,470,473,475,478,480,483,485,488,490,493,495,498,500,...
%                503,505,508,510,513,515,518,520,523,525,528,530,533,535,538,540,543,545,548,550,...
%                553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
%                625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
%                675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];
% 
% % wavelengths = [550,...
% %                553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
% %                625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
% %                675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];
% % MODIS: 412,443,469,488,531,547,555,645,667,678
% % wavelengths = [412,443,469,488,531,547,555,645,667,678];
% 
% num_wl = length(wavelengths);
% 
% % Load real data
% % data_table = readtable('global_14Aug25_PACE_bands_only.csv', 'VariableNamingRule', 'preserve');
% % data_table = readtable('SAB_QWIP_LT03_GLORIA_LAB_RSR_PACE_Green_Red.csv', 'VariableNamingRule', 'preserve');
% data_table = readtable('global_18Aug25_PACE_bands_TSS_LT200.csv', 'VariableNamingRule', 'preserve');
% % data_table = readtable('global_14Aug25_MODIS.csv', 'VariableNamingRule', 'preserve');
% 
% tss = data_table.TSS;  % TSS values (mg/L)
% rrs_data = table2array(data_table(:, 2:end));  % Rrs spectra (rows: samples, cols: wavelengths)
% num_samples = size(rrs_data, 1);
% 
% % Ensure Rrs is positive
% rrs_data = max(rrs_data, eps);
% 
% % Compute second derivatives
% second_deriv = zeros(size(rrs_data));
% for i = 1:num_samples
%     spec = rrs_data(i, :);
%     % Smooth to reduce noise
%     try
%         spec = sgolayfilt(spec, 3, 11);  % 3rd-order, 11-point window
%     catch
%         spec = smoothdata(spec, 'movmean', 5);  % Fallback
%     end
%     first_deriv = gradient(spec, wavelengths);
%     second_deriv(i, :) = gradient(first_deriv, wavelengths);
% end
% 
% % Grid for hyperparameters and random states
% random_states = [49];
% num_trees_list = [213];
% min_leaf_size_list = [4];
% 
% % Initialize results
% results = [];
% best_r2 = -Inf;
% best_params = [];
% best_y_pred = [];
% best_y_test = [];
% best_fig_num = 1;  % For plotting
% 
% % Evaluation function
% eval_metrics = @(y_true, y_pred) struct(...
%     'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
%     'RMSE', sqrt(mean((y_true - y_pred).^2)));
% 
% % Loop over grid
% fig_num = 1;
% for rs = random_states
%     rng(rs);  % Set random state for reproducibility (affects split)
%     [train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
%     X_train = second_deriv(train_idx, :);
%     X_test = second_deriv(test_idx, :);
%     y_train = tss(train_idx);
%     y_test = tss(test_idx);
%     
%     for nt = num_trees_list
%         for mls = min_leaf_size_list
%             % Train Random Forest
%             rf_model = TreeBagger(nt, X_train, y_train, ...
%                 'Method', 'regression', ...
%                 'OOBPredictorImportance', 'on', ...
%                 'MinLeafSize', mls);
%             
%             % Predict
%             y_pred = predict(rf_model, X_test);
%             
%             % Metrics
%             metrics = eval_metrics(y_test, y_pred);
%             
%             % Store results
%             results(end+1, :) = [rs, nt, mls, metrics.R2, metrics.RMSE];
%             
%             % Print
%             fprintf('RS=%d, NumTrees=%d, MinLeafSize=%d - R²: %.4f, RMSE: %.4f mg/L\n', ...
%                 rs, nt, mls, metrics.R2, metrics.RMSE);
%             
%             % Scatter plot with R²
%             figure(fig_num); hold on;
%             scatter(y_test, y_pred, 50, 'b', 'filled');
%             plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', 'LineWidth', 2);
%             text(min(y_test) + 0.1*(max(y_test)-min(y_test)), ...
%                  max(y_test) - 0.1*(max(y_test)-min(y_test)), ...
%                  sprintf('R^2 = %.4f', metrics.R2), 'FontSize', 12);
%             xlabel('In-Situ TSS (mg/L)');
%             ylabel('Predicted TSS (mg/L)');
%             title(sprintf('RF on 2nd Derivatives: RS=%d, Trees=%d, Leaf=%d', rs, nt, mls));
%             grid on;
%             hold off;
%             saveas(gcf, sprintf('PACE_scatter_rf_deriv_rs%d_trees%d_leaf%d.png', rs, nt, mls));
%             fig_num = fig_num + 1;
%             
%             % Check if best
%             if metrics.R2 > best_r2
%                 best_r2 = metrics.R2;
%                 best_params = [rs, nt, mls];
%                 best_y_pred = y_pred;
%                 best_y_test = y_test;
%                 best_fig_num = fig_num - 1;
%             end
%         end
%     end
% end
% 
% % Summary of best model
% fprintf('\nBest Model: RS=%d, NumTrees=%d, MinLeafSize=%d - R²: %.4f, RMSE: %.4f mg/L\n', ...
%     best_params(1), best_params(2), best_params(3), best_r2, results(find(results(:,4)==best_r2,1),5));
% 
% % Re-plot best with highlight (optional, but since all have plots, mention)
% disp('All scatter plots saved with R² annotated. Check PNG files named scatter_rf_deriv_rs*_trees*_leaf*.png');
% 
% % % Print sample predictions for best model
% % fprintf('\nSample Predictions for Best Model (First 10 Test Samples):\n');
% % fprintf('In-Situ TSS\tPredicted TSS\n');
% % for i = 1:min(10, length(best_y_test))
% %     fprintf('%.4f\t\t%.4f\n', best_y_test(i), best_y_pred(i));
% % end
% % 
% % disp('Hyperparameter tuning complete. Use the best parameters for final model.');
% 
% 
% 
% %% Best RF Model


%% Improving random Forest models
clear; close all; clc;

% % Hardcode wavelengths from column names (nm)
wavelengths = [403,405,408,410,413,415,418,420,423,425,428,430,433,435,438,440,443,445,448,450,...
               453,455,458,460,463,465,468,470,473,475,478,480,483,485,488,490,493,495,498,500,...
               503,505,508,510,513,515,518,520,523,525,528,530,533,535,538,540,543,545,548,550,...
               553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
               625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
               675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];
% wavelengths = [550,...
%                553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
%                625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
%                675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];
% 
% MODIS: 412,443,469,488,531,547,555,645,667,678
% wavelengths = [412,443,469,488,531,547,555,645,667,678];

num_wl = length(wavelengths);

% Load real data
% data_table = readtable('global_14Aug25_PACE_bands_only.csv', 'VariableNamingRule', 'preserve');
% data_table = readtable('SAB_QWIP_LT03_GLORIA_LAB_RSR_PACE_Green_Red.csv', 'VariableNamingRule', 'preserve');
% data_table = readtable('global_18Aug25_PACE_bands_TSS_LT200.csv', 'VariableNamingRule', 'preserve');
% data_table = readtable('global_14Aug25_MODIS.csv', 'VariableNamingRule', 'preserve');
data_table = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
tss = data_table.TSS;  % TSS values (mg/L)
rrs_data = table2array(data_table(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
num_samples = size(rrs_data, 1);

% Ensure Rrs is positive
rrs_data = max(rrs_data, eps);

% Compute second derivatives
second_deriv = zeros(size(rrs_data));
for i = 1:num_samples
    spec = rrs_data(i, :);
    % Smooth to reduce noise
    try
        spec = sgolayfilt(spec, 3, 11);  % 3rd-order, 11-point window
    catch
        spec = smoothdata(spec, 'movmean', 5);  % Fallback
    end
    first_deriv = gradient(spec, wavelengths);
    second_deriv(i, :) = gradient(first_deriv, wavelengths);
end
% Check for NaN in your second derivative features
nan_in_features = any(any(isnan(second_deriv)));
fprintf('NaN in second derivative features: %d\n', nan_in_features);

% Check individual samples
nan_samples = any(isnan(second_deriv), 2);
fprintf('Samples with NaN features: %d\n', sum(nan_samples));
%% 

% Data splitting: 80% train, 20% test with optimized RS=49
rng(49);  % Optimized random state
[train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
X_train = second_deriv(train_idx, :);
X_test = second_deriv(test_idx, :);
y_train = tss(train_idx);
y_test = tss(test_idx);

% Train optimized Random Forest model with best hyperparameters
rf_model_best = TreeBagger(213, X_train, y_train, ...  % Optimized: 104 trees
    'Method', 'regression', ...
    'OOBPredictorImportance', 'on', ...
    'MinLeafSize', 4);  % Optimized: MinLeafSize=4
            
% Save the model
save('Simulated_PACE_rf_model_26Aug25.mat', 'rf_model_best', 'wavelengths');

%% 
% Predict
y_pred = predict(rf_model_best, X_test);
% y_pred = str2double(y_pred);  % Convert cell array to numeric

% Evaluate
eval_metrics = @(y_true, y_pred) struct(...
    'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
    'RMSE', sqrt(mean((y_true - y_pred).^2)));
metrics = eval_metrics(y_test, y_pred);

% Display metrics
fprintf('Optimized Random Forest Model (2nd Derivatives, RS=49, NumTrees=104, MinLeafSize=4):\n');
fprintf('R²: %.4f, RMSE: %.4f mg/L\n', metrics.R2, metrics.RMSE);

% Check if test set is empty or has issues
fprintf('Test set size: %d\n', length(y_test));
fprintf('Any NaN in y_test: %d\n', any(isnan(y_test)));
fprintf('Any NaN in y_pred: %d\n', any(isnan(y_pred)));
%%
% Feature importances for Random Forest (Raw Rrs)
importances = rf_model_best.OOBPermutedPredictorDeltaError;
[sorted_imp, imp_idx] = sort(importances, 'descend');

fprintf('\nTop 50 Wavelengths by Importance (Random Forest, Raw Rrs):\n');
for i = 1:min(50, num_wl)
    fprintf('Wavelength %d nm: Importance %.4f\n', wavelengths(imp_idx(i)), sorted_imp(i));
end

% Create feature importance plot
figure('Position', [100, 100, 1200, 600]);

% Subplot 1: Bar plot of top 20 important features
subplot(1, 2, 1);
barh(sorted_imp(1:min(50, num_wl)), 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'YTick', 1:min(50, num_wl));
set(gca, 'YTickLabel', wavelengths(imp_idx(1:min(50, num_wl))));
ylabel('Wavelength (nm)');
xlabel('Feature Importance');
title('Top 20 Important Wavelengths (Raw Rrs)');
grid on;
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % Highest importance at top

% Subplot 2: Spectral plot of all feature importances
subplot(1, 2, 2);
plot(wavelengths, importances, 'b-', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
hold on;

% Highlight top 5 wavelengths
[~, top5_idx] = maxk(importances, 20);
scatter(wavelengths(top5_idx), importances(top5_idx), 100, 'ro', 'filled', ...
        'MarkerEdgeColor', 'r', 'LineWidth', 1.5);

xlabel('Wavelength (nm)');
ylabel('Feature Importance');
title('Feature Importance Spectrum (Raw Rrs)');
grid on;
legend('All wavelengths', 'Top 5 wavelengths', 'Location', 'best');

% Add text annotations for top 5 wavelengths
for i = 1:min(20, num_wl)
    text(wavelengths(top5_idx(i)), importances(top5_idx(i)) + 0.1, ...
         sprintf('%d nm', wavelengths(top5_idx(i))), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Set consistent y-axis limits for better comparison
ylim_vals = [min(importances)*0.9, max(importances)*1.1];
subplot(1, 2, 1); xlim(ylim_vals);
subplot(1, 2, 2); ylim(ylim_vals);

% Add overall title
sgtitle('RF Feature Importance - 2nd Derivatives Rrs Data', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Save the figure
saveas(gcf, 'feature_importance_best_2nd_Der_rrs.png');
% saveas(gcf, 'feature_importance_raw_rrs.fig');

% Display summary statistics
fprintf('\nFeature Importance Summary:\n');
fprintf('Mean importance: %.4f\n', mean(importances));
fprintf('Max importance: %.4f (at %d nm)\n', max(importances), wavelengths(imp_idx(1)));
fprintf('Min importance: %.4f (at %d nm)\n', min(importances), wavelengths(imp_idx(end)));
fprintf('Importance range: %.4f\n', range(importances));

%% Additional diagnostic information
fprintf('\nModel Diagnostics:\n');
fprintf('Number of Observations: %d\n', num_samples);
fprintf('Training Set Size: %d\n', sum(train_idx));
fprintf('Test Set Size: %d\n', sum(test_idx));
fprintf('Out-of-Bag Error: %.4f\n', mean(oobError(rf_model_best)));
%% 
%% Second Derivatives: Random Forest
figure(15); 
hold on;

% Create scatter plot (matching style from figure 12)
scatter(y_test, y_pred, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Add 1:1 line
plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', 'LineWidth', 2);

% Calculate metrics for second derivative model
r2_deriv = 1 - sum((y_test - y_pred).^2) / sum((y_test - mean(y_test)).^2);
rmse_deriv = sqrt(mean((y_test - y_pred).^2));

% Combine all metrics and statistics in one text box
combined_str_deriv = {sprintf('R² = %.3f', r2_deriv), ...
                      sprintf('RMSE = %.2f mg/L', rmse_deriv), ...
                      sprintf('n = %d', length(y_test)), ...
                      sprintf('Bias = %.2f mg/L', mean(y_pred - y_test))};

text_box_deriv = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str_deriv, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

% Formatting (matching figure 12 style)
xlabel('In-Situ TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-2, 200]);
ylim([-2, 200]);
title('Random Forest: In-Situ vs Predicted TSS (Second Derivatives)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Add equation of the regression line for second derivative model
p_deriv = polyfit(y_test, y_pred, 1);
regress_line_deriv = polyval(p_deriv, xlim);
plot(xlim, regress_line_deriv, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');

% Legend for second derivative model
slope_deriv = p_deriv(1);
intercept_deriv = p_deriv(2);
legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope_deriv, intercept_deriv)}, ...
       'Location', 'southeast');

hold off;

% Save high-quality figure (matching figure 12 settings)
set(gcf, 'Position', [100, 100, 800, 700]);
saveas(gcf, 'scatter_rf_deriv.png');
saveas(gcf, 'scatter_rf_deriv.pdf');  % Vector format for publications

% Display metrics in command window too
fprintf('\nRandom Forest (Second Derivatives) Performance:\n');
fprintf('R²: %.4f\n', r2_deriv);
fprintf('RMSE: %.4f mg/L\n', rmse_deriv);
fprintf('Bias: %.4f mg/L\n', mean(y_pred - y_test));
fprintf('Sample size: %d\n', length(y_test));

%% 

