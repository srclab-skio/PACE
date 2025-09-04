
% Load real data
% data_table = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
clear; close all; clc;

% Hardcode wavelengths from column names (nm)
wavelengths = [403,405,408,410,413,415,418,420,423,425,428,430,433,435,438,440,443,445,448,450,...
               453,455,458,460,463,465,468,470,473,475,478,480,483,485,488,490,493,495,498,500,...
               503,505,508,510,513,515,518,520,523,525,528,530,533,535,538,540,543,545,548,550,...
               553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
               625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
               675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];

num_wl = length(wavelengths);

% Load real data
% data_table = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
data_table = readtable('Simul_SAB_Insitu_Rrs_Nechad_et_al_2015_PACE_Rrs_SRS.csv', 'VariableNamingRule', 'preserve');
tss = data_table.TSS;  % TSS values (mg/L)
rrs_data = table2array(data_table(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
num_samples = size(rrs_data, 1);

% Convert response variable to double and check for validity
tss = double(tss);   % Ensure double

% Check for NaN or Inf in response variable
nan_idx = isnan(tss) | isinf(tss);
if any(nan_idx)
    fprintf('Found %d samples with NaN or Inf in TSS. Removing them.\n', sum(nan_idx));
    tss = tss(~nan_idx);
    rrs_data = rrs_data(~nan_idx, :);
    num_samples = size(rrs_data, 1);
end

% Ensure Rrs is positive
rrs_data = max(rrs_data, eps);

% Find index for 665 nm
idx_665 = find(wavelengths == 665);
if isempty(idx_665)
    [~, idx_665] = min(abs(wavelengths - 665));
end

% Compute xtra feature
rrs_665 = rrs_data(:, idx_665);
xtra = 12.74 + (2600 * pi * rrs_665) ./ (1 - (pi * rrs_665) / 1728);
xtra = xtra(:);  % Ensure column vector

% Check for NaN or Inf in xtra
if any(isnan(xtra) | isinf(xtra))
    fprintf('Warning: xtra feature contains NaN or Inf. Replacing with mean.\n');
    xtra(isnan(xtra) | isinf(xtra)) = mean(xtra(~isnan(xtra) & ~isinf(xtra)));
end

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

% Check for NaN in second derivative features
nan_in_features = any(any(isnan(second_deriv)));
fprintf('NaN in second derivative features: %d\n', nan_in_features);

% Check individual samples
nan_samples = any(isnan(second_deriv), 2);
if any(nan_samples)
    fprintf('Found %d samples with NaN in second derivatives. Removing them.\n', sum(nan_samples));
    second_deriv = second_deriv(~nan_samples, :);
    tss = tss(~nan_samples);
    xtra = xtra(~nan_samples);
    num_samples = size(second_deriv, 1);
end

% Hyperparameter grids for tuning
random_states = [ 49];
min_leaf_sizes = [ 4];
num_trees_grid = [ 200];

% Initialize best parameters
best_rmse = inf;
best_rs = [];
best_leaf = [];
best_trees = [];
best_model = [];
best_top20_idx = [];
tic;
% Loop over hyperparameters
for rs = random_states
    rng(rs);  % Set random state for data split
    [train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
    
    X_train_initial = second_deriv(train_idx, :);
    X_test_initial = second_deriv(test_idx, :);
    y_train = tss(train_idx);
    y_test = tss(test_idx);
    
    % Train initial Random Forest model to get feature importances (using default num_trees=100 for speed)
    rf_model_initial = TreeBagger(200, X_train_initial, y_train, ...
        'Method', 'regression', ...
        'OOBPredictorImportance', 'on', ...
        'MinLeafSize', 4);
    
    % Get feature importances from initial model
    importances = rf_model_initial.OOBPermutedPredictorDeltaError;
    [~, imp_idx] = sort(importances, 'descend');
    
    % Select top 20 important features (columns)
    top20_idx = imp_idx(1:20);
    
    % Augment the feature set by duplicating the top 20 features more times (original + 2 duplicates = 3 times total)
    second_deriv_aug = [second_deriv, second_deriv(:, top20_idx), second_deriv(:, top20_idx)];
    
    % Add xtra feature with extra weight (5 times total)
    second_deriv_aug = [second_deriv_aug, repmat(xtra, 1, 5)];
    
    % Update train/test sets with augmented features
    X_train = second_deriv_aug(train_idx, :);
    X_test = second_deriv_aug(test_idx, :);
    
    % Now tune MinLeafSize and NumTrees
    for min_leaf = min_leaf_sizes
        for num_trees = num_trees_grid
            % Train model with current hyperparameters
            rf_model = TreeBagger(num_trees, X_train, y_train, ...
                'Method', 'regression', ...
                'OOBPredictorImportance', 'on', ...
                'MinLeafSize', min_leaf);
            
            % Predict on test set
            y_pred = predict(rf_model, X_test);
            
            % Compute RMSE
            rmse = sqrt(mean((y_test - y_pred).^2));
            
            % Check if this is the best
            if rmse < best_rmse
                best_rmse = rmse;
                best_rs = rs;
                best_leaf = min_leaf;
                best_trees = num_trees;
                best_model = rf_model;
                best_top20_idx = top20_idx;
            end
        end
    end
end

% Display best hyperparameters
fprintf('Best Hyperparameters:\n');
fprintf('Random State: %d\n', best_rs);
fprintf('MinLeafSize: %d\n', best_leaf);
fprintf('NumTrees: %d\n', best_trees);
fprintf('Best RMSE: %.4f mg/L\n', best_rmse);
toc;
% Augment wavelengths for plotting (duplicates labeled with same values, xtra as 0)
wavelengths_aug = [wavelengths, wavelengths(best_top20_idx), wavelengths(best_top20_idx), zeros(1, 5)];

% Save the best weighted model
save('PACE_rf_model_2nd_deriv_weighted_best.mat', 'best_model', 'wavelengths_aug', 'best_top20_idx', 'best_rs', 'best_leaf', 'best_trees');

% Re-split data with best random state for final evaluation and plots
rng(best_rs);
[train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
X_train = second_deriv_aug(train_idx, :);
X_test = second_deriv_aug(test_idx, :);
y_train = tss(train_idx);
y_test = tss(test_idx);

% Since best_model is already trained on this split (from the loop), predict again for confirmation
y_pred = predict(best_model, X_test);

% Evaluate performance
eval_metrics = @(y_true, y_pred) struct(...
    'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
    'RMSE', sqrt(mean((y_true - y_pred).^2)));
metrics = eval_metrics(y_test, y_pred);

% Display metrics
fprintf('Best Optimized Random Forest Model (2nd Derivatives Weighted):\n');
fprintf('R²: %.4f, RMSE: %.4f mg/L\n', metrics.R2, metrics.RMSE);

% Check if test set is empty or has issues
fprintf('Test set size: %d\n', length(y_test));
fprintf('Any NaN in y_test: %d\n', any(isnan(y_test)));
fprintf('Any NaN in y_pred: %d\n', any(isnan(y_pred)));

% Feature importances for best model
importances_aug = best_model.OOBPermutedPredictorDeltaError;
[sorted_imp_aug, imp_idx_aug] = sort(importances_aug, 'descend');

fprintf('\nTop 50 Wavelengths by Importance (Random Forest, 2nd Derivatives Weighted):\n');
for i = 1:min(50, length(wavelengths_aug))
    wl = wavelengths_aug(imp_idx_aug(i));
    if wl == 0
        fprintf('Feature xtra: Importance %.4f\n', sorted_imp_aug(i));
    else
        fprintf('Wavelength %d nm: Importance %.4f\n', wl, sorted_imp_aug(i));
    end
end

% Create feature importance plot (using augmented wavelengths)
figure('Position', [100, 100, 1200, 600]);

% Subplot 1: Bar plot of top 50 important features
subplot(1, 2, 1);
barh(sorted_imp_aug(1:min(50, length(wavelengths_aug))), 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'YTick', 1:min(50, length(wavelengths_aug)));
% Create Y tick labels for top 50 features, labeling 'xtra' for zero wavelengths
ytick_labels = cell(1, min(50, length(wavelengths_aug)));
for i = 1:min(50, length(wavelengths_aug))
    wl = wavelengths_aug(imp_idx_aug(i));
    if wl == 0
        ytick_labels{i} = 'xtra';
    else
        ytick_labels{i} = sprintf('%d nm', wl);
    end
end
set(gca, 'YTickLabel', ytick_labels);
ylabel('Feature');
xlabel('Feature Importance');
title('Top 50 Important Wavelengths (2nd Derivatives Weighted)');
grid on;
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % Highest importance at top

% Subplot 2: Spectral plot of all feature importances
subplot(1, 2, 2);
plot(1:length(wavelengths_aug), importances_aug, 'b-', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);  % Plot by index since duplicates
hold on;

% Highlight top 20 wavelengths (by importance, may include duplicates)
[~, top20_idx_aug] = maxk(importances_aug, 20);
scatter(top20_idx_aug, importances_aug(top20_idx_aug), 100, 'ro', 'filled', ...
        'MarkerEdgeColor', 'r', 'LineWidth', 1.5);

xlabel('Feature Index');
ylabel('Feature Importance');
title('Feature Importance Spectrum (2nd Derivatives Weighted)');
grid on;
legend('All features', 'Top 20 features', 'Location', 'best');

% Add text annotations for top 20 (using wavelengths_aug)
for i = 1:20
    wl = wavelengths_aug(top20_idx_aug(i));
    if wl == 0
        text_str = 'xtra';
    else
        text_str = sprintf('%d nm', wl);
    end
    text(top20_idx_aug(i), importances_aug(top20_idx_aug(i)) + 0.1, ...
         text_str, ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Set consistent y-axis limits for better comparison
ylim_vals = [min(importances_aug)*0.9, max(importances_aug)*1.1];
subplot(1, 2, 1); xlim(ylim_vals);
subplot(1, 2, 2); ylim(ylim_vals);

% Add overall title
sgtitle('RF Feature Importance - 2nd Derivatives Weighted Data', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Save the figure
saveas(gcf, 'feature_importance_2nd_deriv_weighted.png');

% Display summary statistics
fprintf('\nFeature Importance Summary:\n');
fprintf('Mean importance: %.4f\n', mean(importances_aug));
max_wl = wavelengths_aug(imp_idx_aug(1));
if max_wl == 0
    fprintf('Max importance: %.4f (at xtra)\n', max(importances_aug));
else
    fprintf('Max importance: %.4f (at %d nm)\n', max(importances_aug), max_wl);
end
min_wl = wavelengths_aug(imp_idx_aug(end));
if min_wl == 0
    fprintf('Min importance: %.4f (at xtra)\n', min(importances_aug));
else
    fprintf('Min importance: %.4f (at %d nm)\n', min(importances_aug), min_wl);
end
fprintf('Importance range: %.4f\n', range(importances_aug));

% Additional diagnostic information
fprintf('\nModel Diagnostics:\n');
fprintf('Number of Observations: %d\n', num_samples);
fprintf('Training Set Size: %d\n', sum(train_idx));
fprintf('Test Set Size: %d\n', sum(test_idx));
fprintf('Out-of-Bag Error: %.4f\n', mean(oobError(best_model)));

% Scatter plot for weighted 2nd derivatives model
figure(15); 
hold on;

% Create scatter plot
scatter(y_test, y_pred, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Add 1:1 line
plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', 'LineWidth', 2);

% Calculate metrics
r2_weighted = 1 - sum((y_test - y_pred).^2) / sum((y_test - mean(y_test)).^2);
rmse_weighted = sqrt(mean((y_test - y_pred).^2));

% Combine metrics and statistics in one text box
combined_str_weighted = {sprintf('R² = %.3f', r2_weighted), ...
                         sprintf('RMSE = %.2f mg/L', rmse_weighted), ...
                         sprintf('n = %d', length(y_test)), ...
                         sprintf('Bias = %.2f mg/L', mean(y_pred - y_test))};

text_box_weighted = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str_weighted, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

% Formatting
xlabel('In-Situ TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-2, 200]);
ylim([-2, 200]);
title('Random Forest: In-Situ vs Predicted TSS (2nd Derivatives Weighted)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Add equation of the regression line
p_weighted = polyfit(y_test, y_pred, 1);
regress_line_weighted = polyval(p_weighted, xlim);
plot(xlim, regress_line_weighted, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');

% Legend
slope_weighted = p_weighted(1);
intercept_weighted = p_weighted(2);
legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope_weighted, intercept_weighted)}, ...
       'Location', 'southeast');

hold off;

% Save high-quality figure
set(gcf, 'Position', [100, 100, 800, 700]);
saveas(gcf, 'scatter_rf_2nd_deriv_weighted.png');

% Display metrics in command window
fprintf('\nRandom Forest (2nd Derivatives Weighted) Performance:\n');
fprintf('R²: %.4f\n', r2_weighted);
fprintf('RMSE: %.4f mg/L\n', rmse_weighted);
fprintf('Bias: %.4f mg/L\n', mean(y_pred - y_test));
fprintf('Sample size: %d\n', length(y_test));

%% New section
%% The best is RS= 42, ntrees = 213, nleafs = 4
clear; close all; clc;

% Hardcode wavelengths from column names (nm)
wavelengths = [403,405,408,410,413,415,418,420,423,425,428,430,433,435,438,440,443,445,448,450,...
               453,455,458,460,463,465,468,470,473,475,478,480,483,485,488,490,493,495,498,500,...
               503,505,508,510,513,515,518,520,523,525,528,530,533,535,538,540,543,545,548,550,...
               553,555,558,560,563,565,568,570,573,575,578,580,583,585,588,613,615,618,620,623,...
               625,628,630,633,635,638,640,643,645,648,650,653,655,658,660,663,665,668,670,673,...
               675,678,680,683,685,688,690,693,695,698,700,703,705,708,710,713,718];

num_wl = length(wavelengths);

% Load real data
data_table = readtable('Simul_SAB_Insitu_Rrs_Nechad_et_al_2015_PACE_Rrs_SRS.csv', 'VariableNamingRule', 'preserve');
tss = data_table.TSS;  % TSS values (mg/L)
rrs_data = table2array(data_table(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
num_samples = size(rrs_data, 1);

% Convert response variable to double and check for validity
tss = double(tss);   % Ensure double

% Check for NaN or Inf in response variable
nan_idx = isnan(tss) | isinf(tss);
if any(nan_idx)
    fprintf('Found %d samples with NaN or Inf in TSS. Removing them.\n', sum(nan_idx));
    tss = tss(~nan_idx);
    rrs_data = rrs_data(~nan_idx, :);
    num_samples = size(rrs_data, 1);
end

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

% Check for NaN in second derivative features
nan_in_features = any(any(isnan(second_deriv)));
fprintf('NaN in second derivative features: %d\n', nan_in_features);

% Check individual samples
nan_samples = any(isnan(second_deriv), 2);
if any(nan_samples)
    fprintf('Found %d samples with NaN in second derivatives. Removing them.\n', sum(nan_samples));
    second_deriv = second_deriv(~nan_samples, :);
    tss = tss(~nan_samples);
    num_samples = size(second_deriv, 1);
end

% Hyperparameter grids for tuning
random_states = [42];
min_leaf_sizes = [4];
num_trees_grid = [200,213];

% Initialize best parameters
best_rmse = inf;
best_rs = [];
best_leaf = [];
best_trees = [];
best_model = [];

% Loop over hyperparameters
for rs = random_states
    rng(rs);  % Set random state for data split
    [train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
    
    X_train = second_deriv(train_idx, :);
    X_test = second_deriv(test_idx, :);
    y_train = tss(train_idx);
    y_test = tss(test_idx);
    
    % Tune MinLeafSize and NumTrees
    for min_leaf = min_leaf_sizes
        for num_trees = num_trees_grid
            % Train model with current hyperparameters
            rf_model = TreeBagger(num_trees, X_train, y_train, ...
                'Method', 'regression', ...
                'OOBPredictorImportance', 'on', ...
                'MinLeafSize', min_leaf);
            
            % Predict on test set
            y_pred = predict(rf_model, X_test);
            
            % Compute RMSE
            rmse = sqrt(mean((y_test - y_pred).^2));
            
            % Check if this is the best
            if rmse < best_rmse
                best_rmse = rmse;
                best_rs = rs;
                best_leaf = min_leaf;
                best_trees = num_trees;
                best_model = rf_model;
            end
        end
    end
end

% Display best hyperparameters
fprintf('Best Hyperparameters:\n');
fprintf('Random State: %d\n', best_rs);
fprintf('MinLeafSize: %d\n', best_leaf);
fprintf('NumTrees: %d\n', best_trees);
fprintf('Best RMSE: %.4f mg/L\n', best_rmse);

% Save the best model
save('PACE_rf_model_2nd_deriv_best_Simulated_insitu.mat', 'best_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');

%% Re-split data with best random state for final evaluation and plots
rng(best_rs);
[train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
X_train = second_deriv(train_idx, :);
X_test = second_deriv(test_idx, :);
y_train = tss(train_idx);
y_test = tss(test_idx);

% Since best_model is already trained on this split (from the loop), predict again for confirmation
y_pred = predict(best_model, X_test);

% Evaluate performance
eval_metrics = @(y_true, y_pred) struct(...
    'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
    'RMSE', sqrt(mean((y_true - y_pred).^2)));
metrics = eval_metrics(y_test, y_pred);

% Display metrics
fprintf('Best Optimized Random Forest Model (2nd Derivatives):\n');
fprintf('R²: %.4f, RMSE: %.4f mg/L\n', metrics.R2, metrics.RMSE);

% Check if test set is empty or has issues
fprintf('Test set size: %d\n', length(y_test));
fprintf('Any NaN in y_test: %d\n', any(isnan(y_test)));
fprintf('Any NaN in y_pred: %d\n', any(isnan(y_pred)));

% Feature importances for best model
importances = best_model.OOBPermutedPredictorDeltaError;
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
set(gca, 'YDir', 'reverse');  % Highest importance at top

% Subplot 2: Spectral plot of all feature importances
subplot(1, 2, 2);
plot(wavelengths, importances, 'b-', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
hold on;

% Highlight top 20 wavelengths
[~, top20_idx] = maxk(importances, 20);
scatter(wavelengths(top20_idx), importances(top20_idx), 100, 'ro', 'filled', ...
        'MarkerEdgeColor', 'r', 'LineWidth', 1.5);

xlabel('Wavelength (nm)');
ylabel('Feature Importance');
title('Feature Importance Spectrum (2nd Derivatives)');
grid on;
legend('All features', 'Top 20 features', 'Location', 'best');

% Add text annotations for top 20
for i = 1:20
    text(wavelengths(top20_idx(i)), importances(top20_idx(i)) + 0.1, ...
         sprintf('%d nm', wavelengths(top20_idx(i))), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Set consistent y-axis limits
ylim_vals = [min(importances)*0.9, max(importances)*1.1];
subplot(1, 2, 1); xlim(ylim_vals);
subplot(1, 2, 2); ylim(ylim_vals);

% Add overall title
sgtitle('RF Feature Importance - 2nd Derivatives Data', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Save the figure
saveas(gcf, 'feature_importance_2nd_deriv.png');

% Display summary statistics
fprintf('\nFeature Importance Summary:\n');
fprintf('Mean importance: %.4f\n', mean(importances));
fprintf('Max importance: %.4f (at %d nm)\n', max(importances), wavelengths(imp_idx(1)));
fprintf('Min importance: %.4f (at %d nm)\n', min(importances), wavelengths(imp_idx(end)));
fprintf('Importance range: %.4f\n', range(importances));

% Additional diagnostic information
fprintf('\nModel Diagnostics:\n');
fprintf('Number of Observations: %d\n', num_samples);
fprintf('Training Set Size: %d\n', sum(train_idx));
fprintf('Test Set Size: %d\n', sum(test_idx));
fprintf('Out-of-Bag Error: %.4f\n', mean(oobError(best_model)));

% Scatter plot for 2nd derivatives model
figure(15); 
hold on;

% Create scatter plot
scatter(y_test, y_pred, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Add 1:1 line
plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', 'LineWidth', 2);

% Calculate metrics
r2 = 1 - sum((y_test - y_pred).^2) / sum((y_test - mean(y_test)).^2);
rmse = sqrt(mean((y_test - y_pred).^2));

% Combine metrics and statistics in one text box
combined_str = {sprintf('R² = %.3f', r2), ...
                sprintf('RMSE = %.2f mg/L', rmse), ...
                sprintf('n = %d', length(y_test)), ...
                sprintf('Bias = %.2f mg/L', mean(y_pred - y_test))};

text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

% Formatting
xlabel('In-Situ TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-2, 200]);
ylim([-2, 200]);
title('Random Forest: In-Situ vs Predicted TSS (2nd Derivatives)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Add equation of the regression line
p = polyfit(y_test, y_pred, 1);
regress_line = polyval(p, xlim);
plot(xlim, regress_line, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');

% Legend
slope = p(1);
intercept = p(2);
legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope, intercept)}, ...
       'Location', 'southeast');

hold off;

% Save high-quality figure
set(gcf, 'Position', [100, 100, 800, 700]);
saveas(gcf, 'scatter_rf_2nd_deriv.png');

% Display metrics in command window
fprintf('\nRandom Forest (2nd Derivatives) Performance:\n');
fprintf('R²: %.4f\n', r2);
fprintf('RMSE: %.4f mg/L\n', rmse);
fprintf('Bias: %.4f mg/L\n', mean(y_pred - y_test));
fprintf('Sample size: %d\n', length(y_test));


%% Increase ranges with simulated data
%% Random Forest with Square Root Transformation and NIR Band Weighting
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
% num_wl = length(wavelengths);
% 
% % Load real data
% data_table = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
% tss = data_table.TSS;  % TSS values (mg/L)
% rrs_data = table2array(data_table(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
% num_samples = size(rrs_data, 1);
% 
% % Convert response variable to double and check for validity
% tss = double(tss);   % Ensure double
% 
% % Ensure TSS is non-negative for square root
% if any(tss <= 0)
%     fprintf('Found %d non-positive TSS values. Setting them to 1e-6.\n', sum(tss <= 0));
%     tss(tss <= 0) = 1e-6;
% end
% 
% % Check for NaN or Inf in response variable
% nan_idx = isnan(tss) | isinf(tss);
% if any(nan_idx)
%     fprintf('Found %d samples with NaN or Inf in TSS. Removing them.\n', sum(nan_idx));
%     tss = tss(~nan_idx);
%     rrs_data = rrs_data(~nan_idx, :);
%     num_samples = size(rrs_data, 1);
% end
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
% % Check for NaN in second derivative features
% nan_in_features = any(any(isnan(second_deriv)));
% fprintf('NaN in second derivative features: %d\n', nan_in_features);
% 
% % Check individual samples
% nan_samples = any(isnan(second_deriv), 2);
% if any(nan_samples)
%     fprintf('Found %d samples with NaN in second derivatives. Removing them.\n', sum(nan_samples));
%     second_deriv = second_deriv(~nan_samples, :);
%     tss = tss(~nan_samples);
%     num_samples = size(second_deriv, 1);
% end
% 
% % Apply NIR weighting (700-718 nm)
% nir_idx = wavelengths >= 700;  % NIR bands
% nir_weight = 2;  % Weighting factor for NIR bands
% fprintf('Applying weight %.2f to %d NIR bands (700-718 nm)\n', nir_weight, sum(nir_idx));
% second_deriv(:, nir_idx) = second_deriv(:, nir_idx) * nir_weight;
% 
% % Hyperparameter grids for tuning
% random_states = [42];
% min_leaf_sizes = [4];
% num_trees_grid = [100, 200, 300];
% 
% % Initialize best parameters
% best_rmse = inf;
% best_rs = [];
% best_leaf = [];
% best_trees = [];
% best_model = [];
% 
% % Loop over hyperparameters
% for rs = random_states
%     rng(rs);  % Set random state for data split
%     [train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
%     
%     X_train = second_deriv(train_idx, :);
%     X_test = second_deriv(test_idx, :);
%     y_train = tss(train_idx);
%     y_test = tss(test_idx);
%     
%     % Apply square root transformation
%     y_train_trans = sqrt(y_train);
%     
%     % Oversample high TSS samples (above 75th percentile, 1x duplication)
%     high_tss_threshold = prctile(y_train, 75);
%     high_tss_idx = y_train > high_tss_threshold;
%     num_high_tss = sum(high_tss_idx);
%     fprintf('Oversampling %d samples with TSS > %.2f mg/L (75th percentile)\n', num_high_tss, high_tss_threshold);
%     
%     % Duplicate high TSS samples (1x)
%     X_train_oversampled = [X_train; X_train(high_tss_idx, :)];
%     y_train_trans_oversampled = [y_train_trans; y_train_trans(high_tss_idx)];
%     
%     % Tune MinLeafSize and NumTrees
%     for min_leaf = min_leaf_sizes
%         for num_trees = num_trees_grid
%             % Train random forest model on square root transformed, oversampled data
%             rf_model = TreeBagger(num_trees, X_train_oversampled, y_train_trans_oversampled, ...
%                 'Method', 'regression', ...
%                 'OOBPrediction', 'on', ...
%                 'MinLeafSize', min_leaf);
%             
%             % Predict on test set (transformed scale)
%             y_pred_trans = predict(rf_model, X_test);
%             
%             % Inverse square root transformation
%             y_pred = y_pred_trans .^ 2;
%             
%             % Compute RMSE on original scale
%             rmse = sqrt(mean((y_test - y_pred).^2));
%             
%             % Check if this is the best
%             if rmse < best_rmse
%                 best_rmse = rmse;
%                 best_rs = rs;
%                 best_leaf = min_leaf;
%                 best_trees = num_trees;
%                 best_model = rf_model;
%             end
%         end
%     end
% end
% 
% % Display best hyperparameters
% fprintf('Best Hyperparameters:\n');
% fprintf('Random State: %d\n', best_rs);
% fprintf('MinLeafSize: %d\n', best_leaf);
% fprintf('NumTrees: %d\n', best_trees);
% fprintf('Best RMSE: %.4f mg/L\n', best_rmse);
% 
% % Save the best model
% save('PACE_rf_model_2nd_deriv_best_Simulated_insitu_sqrt_nirweighted.mat', 'best_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees', 'nir_idx', 'nir_weight');
% 
% %% Re-split data with best random state for final evaluation and plots
% rng(best_rs);
% [train_idx, test_idx] = crossvalind('HoldOut', num_samples, 0.2);
% X_train = second_deriv(train_idx, :);
% X_test = second_deriv(test_idx, :);
% y_train = tss(train_idx);
% y_test = tss(test_idx);
% 
% % Apply square root transformation and oversampling for confirmation
% y_train_trans = sqrt(y_train);
% high_tss_threshold = prctile(y_train, 75);
% high_tss_idx = y_train > high_tss_threshold;
% X_train_oversampled = [X_train; X_train(high_tss_idx, :)];
% y_train_trans_oversampled = [y_train_trans; y_train_trans(high_tss_idx)];
% 
% % Since best_model is already trained, predict again for confirmation
% y_pred_trans = predict(best_model, X_test);
% y_pred = y_pred_trans .^ 2;
% 
% % Evaluate performance
% eval_metrics = @(y_true, y_pred) struct(...
%     'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
%     'RMSE', sqrt(mean((y_true - y_pred).^2)));
% metrics = eval_metrics(y_test, y_pred);
% 
% % Display metrics
% fprintf('Best Optimized Random Forest Model (2nd Derivatives, Sqrt + NIR Weighting):\n');
% fprintf('R²: %.4f, RMSE: %.4f mg/L\n', metrics.R2, metrics.RMSE);
% 
% % Check if test set is empty or has issues
% fprintf('Test set size: %d\n', length(y_test));
% fprintf('Any NaN in y_test: %d\n', any(isnan(y_test)));
% fprintf('Any NaN in y_pred: %d\n', any(isnan(y_pred)));
% 
% % Manual feature importance via permutation
% rng('default');  % For reproducibility
% importances = zeros(1, num_wl);
% y_pred_base = predict(best_model, X_test);
% mse_base = mean((y_test - y_pred_base .^ 2).^2);
% for i = 1:num_wl
%     X_test_perm = X_test;
%     X_test_perm(:, i) = X_test_perm(randperm(size(X_test, 1)), i);  % Shuffle feature i
%     y_pred_perm = predict(best_model, X_test_perm) .^ 2;
%     mse_perm = mean((y_test - y_pred_perm).^2);
%     importances(i) = mse_perm - mse_base;  % Increase in MSE
% end
% [sorted_imp, imp_idx] = sort(importances, 'descend');
% 
% fprintf('\nTop 50 Wavelengths by Importance (Random Forest, 2nd Derivatives, Sqrt + NIR Weighting):\n');
% for i = 1:min(50, length(wavelengths))
%     fprintf('Wavelength %d nm: Importance %.4f\n', wavelengths(imp_idx(i)), sorted_imp(i));
% end
% 
% % Create feature importance plot
% figure('Position', [100, 100, 1200, 600]);
% 
% % Subplot 1: Bar plot of top 50 important features
% subplot(1, 2, 1);
% barh(sorted_imp(1:min(50, length(wavelengths))), 'FaceColor', [0.2, 0.6, 0.8]);
% set(gca, 'YTick', 1:min(50, length(wavelengths)));
% set(gca, 'YTickLabel', wavelengths(imp_idx(1:min(50, length(wavelengths)))));
% ylabel('Wavelength (nm)');
% xlabel('Feature Importance (MSE Increase)');
% title('Top 50 Important Wavelengths (2nd Derivatives, Sqrt + NIR Weighting)');
% grid on;
% set(gca, 'XAxisLocation', 'top');
% set(gca, 'YDir', 'reverse');  % Highest importance at top
% 
% % Subplot 2: Spectral plot of all feature importances
% subplot(1, 2, 2);
% plot(wavelengths, importances, 'b-', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
% hold on;
% 
% % Highlight top 20 wavelengths
% [~, top20_idx] = maxk(importances, 20);
% scatter(wavelengths(top20_idx), importances(top20_idx), 100, 'ro', 'filled', ...
%         'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
% 
% xlabel('Wavelength (nm)');
% ylabel('Feature Importance (MSE Increase)');
% title('Feature Importance Spectrum (2nd Derivatives, Sqrt + NIR Weighting)');
% grid on;
% legend('All features', 'Top 20 features', 'Location', 'best');
% 
% % Add text annotations for top 20
% for i = 1:20
%     text(wavelengths(top20_idx(i)), importances(top20_idx(i)) + 0.1, ...
%          sprintf('%d nm', wavelengths(top20_idx(i))), ...
%          'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% end
% 
% % Set consistent y-axis limits
% ylim_vals = [min(importances)*0.9, max(importances)*1.1];
% subplot(1, 2, 1); xlim(ylim_vals);
% subplot(1, 2, 2); ylim(ylim_vals);
% 
% % Add overall title
% sgtitle('Random Forest Feature Importance - 2nd Derivatives Data (Sqrt + NIR Weighting)', ...
%         'FontSize', 14, 'FontWeight', 'bold');
% 
% % Save the figure
% saveas(gcf, 'feature_importance_2nd_deriv_sqrt_nirweighted.png');
% 
% % Display summary statistics
% fprintf('\nFeature Importance Summary:\n');
% fprintf('Mean importance: %.4f\n', mean(importances));
% fprintf('Max importance: %.4f (at %d nm)\n', max(importances), wavelengths(imp_idx(1)));
% fprintf('Min importance: %.4f (at %d nm)\n', min(importances), wavelengths(imp_idx(end)));
% fprintf('Importance range: %.4f\n', range(importances));
% 
% % Additional diagnostic information
% fprintf('\nModel Diagnostics:\n');
% fprintf('Number of Observations: %d\n', num_samples);
% fprintf('Training Set Size (after oversampling): %d\n', size(X_train_oversampled, 1));
% fprintf('Test Set Size: %d\n', sum(test_idx));
% fprintf('Out-of-Bag Error: %.4f\n', mean(oobError(best_model)));
% 
% % Scatter plot for predictions
% figure(15); 
% hold on;
% 
% % Create scatter plot
% scatter(y_test, y_pred, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
% 
% % Add 1:1 line
% plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', 'LineWidth', 2);
% 
% % Calculate metrics
% r2 = 1 - sum((y_test - y_pred).^2) / sum((y_test - mean(y_test)).^2);
% rmse = sqrt(mean((y_test - y_pred).^2));
% 
% % Combine metrics and statistics in one text box
% combined_str = {sprintf('R² = %.3f', r2), ...
%                 sprintf('RMSE = %.2f mg/L', rmse), ...
%                 sprintf('n = %d', length(y_test)), ...
%                 sprintf('Bias = %.2f mg/L', mean(y_pred - y_test))};
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
% title('Random Forest: In-Situ vs Predicted TSS (2nd Derivatives, Sqrt + NIR Weighting)', 'FontSize', 14, 'FontWeight', 'bold');
% grid on;
% set(gca, 'FontSize', 11);
% 
% % Add equation of the regression line
% p = polyfit(y_test, y_pred, 1);
% regress_line = polyval(p, xlim);
% plot(xlim, regress_line, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');
% 
% % Legend
% slope = p(1);
% intercept = p(2);
% legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope, intercept)}, ...
%        'Location', 'southeast');
% 
% hold off;
% 
% % Save high-quality figure
% set(gcf, 'Position', [100, 100, 800, 700]);
% saveas(gcf, 'scatter_rf_2nd_deriv_sqrt_nirweighted.png');
% 
% % Display metrics in command window
% fprintf('\nRandom Forest (2nd Derivatives, Sqrt + NIR Weighting) Performance:\n');
% fprintf('R²: %.4f\n', r2);
% fprintf('RMSE: %.4f mg/L\n', rmse);
% fprintf('Bias: %.4f mg/L\n', mean(y_pred - y_test));
% fprintf('Sample size: %d\n', length(y_test));


%% train with simulated and test with insitu Gloria +Lab
% %% Random Forest with Second Derivatives, Validated on In-Situ Data
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
% num_wl = length(wavelengths);
% 
% % Load simulated data (training)
% data_table_sim = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
% tss_sim = data_table_sim.TSS;  % TSS values (mg/L)
% rrs_data_sim = table2array(data_table_sim(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
% num_samples_sim = size(rrs_data_sim, 1);
% 
% % Convert response variable to double and check for validity
% tss_sim = double(tss_sim);
% % Check for NaN or Inf in response variable
% nan_idx_sim = isnan(tss_sim) | isinf(tss_sim);
% if any(nan_idx_sim)
%     fprintf('Found %d samples with NaN or Inf in TSS (simulated). Removing them.\n', sum(nan_idx_sim));
%     tss_sim = tss_sim(~nan_idx_sim);
%     rrs_data_sim = rrs_data_sim(~nan_idx_sim, :);
%     num_samples_sim = size(rrs_data_sim, 1);
% end
% % Ensure Rrs is positive
% rrs_data_sim = max(rrs_data_sim, eps);
% 
% % Load in-situ data (testing, rows 465:end)
% data_table_insitu = readtable('global_14Aug25_PACE_bands_only.csv', 'VariableNamingRule', 'preserve');
% tss_insitu = data_table_insitu.TSS;  % TSS values (mg/L)
% rrs_data_insitu = table2array(data_table_insitu(:, 4:end));  % Rrs spectra
% num_samples_insitu = size(rrs_data_insitu, 1);
% 
% % Convert response variable to double and check for validity
% tss_insitu = double(tss_insitu);
% % Check for NaN or Inf in response variable
% nan_idx_insitu = isnan(tss_insitu) | isinf(tss_insitu);
% if any(nan_idx_insitu)
%     fprintf('Found %d samples with NaN or Inf in TSS (in-situ). Removing them.\n', sum(nan_idx_insitu));
%     tss_insitu = tss_insitu(~nan_idx_insitu);
%     rrs_data_insitu = rrs_data_insitu(~nan_idx_insitu, :);
%     num_samples_insitu = size(rrs_data_insitu, 1);
% end
% % Ensure Rrs is positive
% rrs_data_insitu = max(rrs_data_insitu, eps);
% 
% % Interpolate in-situ data to model wavelengths (assuming in-situ wavelengths are 400-800 nm)
% insitu_wl = 400:800;  % Placeholder; replace with actual wavelengths if available
% if size(rrs_data_insitu, 2) ~= num_wl
%     fprintf('In-situ data has %d wavelengths; interpolating to match model wavelengths (403-718 nm).\n', size(rrs_data_insitu, 2));
%     rrs_data_insitu_interp = zeros(num_samples_insitu, num_wl);
%     for i = 1:num_samples_insitu
%         rrs_data_insitu_interp(i, :) = interp1(insitu_wl, rrs_data_insitu(i, :), wavelengths, 'linear', 'extrap');
%     end
%     rrs_data_insitu = rrs_data_insitu_interp;
% end
% 
% % Compute second derivatives for simulated data (training)
% second_deriv_sim = zeros(size(rrs_data_sim));
% for i = 1:num_samples_sim
%     spec = rrs_data_sim(i, :);
%     try
%         spec = sgolayfilt(spec, 3, 11);  % 3rd-order, 11-point window
%     catch
%         spec = smoothdata(spec, 'movmean', 5);  % Fallback
%     end
%     first_deriv = gradient(spec, wavelengths);
%     second_deriv_sim(i, :) = gradient(first_deriv, wavelengths);
% end
% 
% % Check for NaN in second derivative features (simulated)
% nan_in_features_sim = any(any(isnan(second_deriv_sim)));
% fprintf('NaN in second derivative features (simulated): %d\n', nan_in_features_sim);
% nan_samples_sim = any(isnan(second_deriv_sim), 2);
% if any(nan_samples_sim)
%     fprintf('Found %d samples with NaN in second derivatives (simulated). Removing them.\n', sum(nan_samples_sim));
%     second_deriv_sim = second_deriv_sim(~nan_samples_sim, :);
%     tss_sim = tss_sim(~nan_samples_sim);
%     num_samples_sim = size(second_deriv_sim, 1);
% end
% 
% % Compute second derivatives for in-situ data (testing)
% second_deriv_insitu = zeros(size(rrs_data_insitu));
% for i = 1:num_samples_insitu
%     spec = rrs_data_insitu(i, :);
%     try
%         spec = sgolayfilt(spec, 3, 11);
%     catch
%         spec = smoothdata(spec, 'movmean', 5);
%     end
%     first_deriv = gradient(spec, wavelengths);
%     second_deriv_insitu(i, :) = gradient(first_deriv, wavelengths);
% end
% 
% % Check for NaN in second derivative features (in-situ)
% nan_in_features_insitu = any(any(isnan(second_deriv_insitu)));
% fprintf('NaN in second derivative features (in-situ): %d\n', nan_in_features_insitu);
% nan_samples_insitu = any(isnan(second_deriv_insitu), 2);
% if any(nan_samples_insitu)
%     fprintf('Found %d samples with NaN in second derivatives (in-situ). Removing them.\n', sum(nan_samples_insitu));
%     second_deriv_insitu = second_deriv_insitu(~nan_samples_insitu, :);
%     tss_insitu = tss_insitu(~nan_samples_insitu);
%     num_samples_insitu = size(second_deriv_insitu, 1);
% end
% 
% % Use best hyperparameters from original code
% best_rs = 42;
% best_leaf = 4;
% best_trees = 213;
% 
% % Train model on full simulated data
% rng(best_rs);
% rf_model = TreeBagger(best_trees, second_deriv_sim, tss_sim, ...
%     'Method', 'regression', ...
%     'MinLeafSize', best_leaf);
% 
% % Save the model
% save('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', ...
%      'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% 
% % Predict on in-situ test data
% y_pred = predict(rf_model, second_deriv_insitu);
% 
% % Evaluate performance
% eval_metrics = @(y_true, y_pred) struct(...
%     'R2', 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2), ...
%     'RMSE', sqrt(mean((y_true - y_pred).^2)));
% metrics = eval_metrics(tss_insitu, y_pred);
% 
% % Additional metrics for high TSS (>50 mg/L)
% high_tss_idx = tss_insitu > 50;
% metrics_high = eval_metrics(tss_insitu(high_tss_idx), y_pred(high_tss_idx));
% 
% % Display metrics
% fprintf('Random Forest Model (2nd Derivatives, In-Situ Validation):\n');
% fprintf('Overall R²: %.4f, RMSE: %.4f mg/L\n', metrics.R2, metrics.RMSE);
% fprintf('High TSS (>50 mg/L) R²: %.4f, RMSE: %.4f mg/L, n: %d\n', ...
%         metrics_high.R2, metrics_high.RMSE, sum(high_tss_idx));
% 
% % Check test set diagnostics
% fprintf('Test set size (in-situ): %d\n', length(tss_insitu));
% fprintf('Any NaN in y_test: %d\n', any(isnan(tss_insitu)));
% fprintf('Any NaN in y_pred: %d\n', any(isnan(y_pred)));
% 
% % Manual feature importance via permutation
% rng('default');
% importances = zeros(1, num_wl);
% y_pred_base = predict(rf_model, second_deriv_insitu);
% mse_base = mean((tss_insitu - y_pred_base).^2);
% for i = 1:num_wl
%     X_test_perm = second_deriv_insitu;
%     X_test_perm(:, i) = X_test_perm(randperm(size(second_deriv_insitu, 1)), i);
%     y_pred_perm = predict(rf_model, X_test_perm);
%     mse_perm = mean((tss_insitu - y_pred_perm).^2);
%     importances(i) = mse_perm - mse_base;  % Increase in MSE
% end
% [sorted_imp, imp_idx] = sort(importances, 'descend');
% 
% fprintf('\nTop 50 Wavelengths by Importance (Random Forest, 2nd Derivatives):\n');
% for i = 1:min(50, length(wavelengths))
%     fprintf('Wavelength %d nm: Importance %.4f\n', wavelengths(imp_idx(i)), sorted_imp(i));
% end
% 
% % Create feature importance plot
% figure('Position', [100, 100, 1200, 600]);
% 
% % Subplot 1: Bar plot of top 50 important features
% subplot(1, 2, 1);
% barh(sorted_imp(1:min(50, length(wavelengths))), 'FaceColor', [0.2, 0.6, 0.8]);
% set(gca, 'YTick', 1:min(50, length(wavelengths)));
% set(gca, 'YTickLabel', wavelengths(imp_idx(1:min(50, length(wavelengths)))));
% ylabel('Wavelength (nm)');
% xlabel('Feature Importance (MSE Increase)');
% title('Top 50 Important Wavelengths (2nd Derivatives)');
% grid on;
% set(gca, 'XAxisLocation', 'top');
% set(gca, 'YDir', 'reverse');
% 
% % Subplot 2: Spectral plot of all feature importances
% subplot(1, 2, 2);
% plot(wavelengths, importances, 'b-', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
% hold on;
% [~, top20_idx] = maxk(importances, 20);
% scatter(wavelengths(top20_idx), importances(top20_idx), 100, 'ro', 'filled', ...
%         'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
% xlabel('Wavelength (nm)');
% ylabel('Feature Importance (MSE Increase)');
% title('Feature Importance Spectrum (2nd Derivatives)');
% grid on;
% legend('All features', 'Top 20 features', 'Location', 'best');
% for i = 1:20
%     text(wavelengths(top20_idx(i)), importances(top20_idx(i)) + 0.1, ...
%          sprintf('%d nm', wavelengths(top20_idx(i))), ...
%          'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% end
% ylim_vals = [min(importances)*0.9, max(importances)*1.1];
% subplot(1, 2, 1); xlim(ylim_vals);
% subplot(1, 2, 2); ylim(ylim_vals);
% sgtitle('Random Forest Feature Importance - 2nd Derivatives (In-Situ Validation)', ...
%         'FontSize', 14, 'FontWeight', 'bold');
% saveas(gcf, 'feature_importance_2nd_deriv_validated.png');
% 
% % Display feature importance summary
% fprintf('\nFeature Importance Summary:\n');
% fprintf('Mean importance: %.4f\n', mean(importances));
% fprintf('Max importance: %.4f (at %d nm)\n', max(importances), wavelengths(imp_idx(1)));
% fprintf('Min importance: %.4f (at %d nm)\n', min(importances), wavelengths(imp_idx(end)));
% fprintf('Importance range: %.4f\n', range(importances));
% 
% % Model diagnostics
% fprintf('\nModel Diagnostics:\n');
% fprintf('Number of Observations (simulated, training): %d\n', num_samples_sim);
% fprintf('Number of Observations (in-situ, testing): %d\n', num_samples_insitu);
% 
% % Scatter plot for predictions
% figure(15);
% hold on;
% scatter(tss_insitu, y_pred, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
% plot([min(tss_insitu), max(tss_insitu)], [min(tss_insitu), max(tss_insitu)], 'r--', 'LineWidth', 2);
% r2 = 1 - sum((tss_insitu - y_pred).^2) / sum((tss_insitu - mean(tss_insitu)).^2);
% rmse = sqrt(mean((tss_insitu - y_pred).^2));
% combined_str = {sprintf('R² = %.3f', r2), ...
%                 sprintf('RMSE = %.2f mg/L', rmse), ...
%                 sprintf('n = %d', length(tss_insitu)), ...
%                 sprintf('Bias = %.2f mg/L', mean(y_pred - tss_insitu))};
% text_box = annotation('textbox', [0.15, 0.75, 0.2, 0.12], 'String', combined_str, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
%     'FontSize', 12, 'FontWeight', 'bold', 'Margin', 5, ...
%     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
% xlabel('In-Situ TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Predicted TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
% xlim([-2, 200]);
% ylim([-2, 200]);
% title('Random Forest: In-Situ vs Predicted TSS (2nd Derivatives, In-Situ Validation)', ...
%       'FontSize', 14, 'FontWeight', 'bold');
% grid on;
% set(gca, 'FontSize', 11);
% p = polyfit(tss_insitu, y_pred, 1);
% regress_line = polyval(p, xlim);
% plot(xlim, regress_line, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');
% slope = p(1);
% intercept = p(2);
% legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope, intercept)}, ...
%        'Location', 'southeast');
% hold off;
% set(gcf, 'Position', [100, 100, 800, 700]);
% saveas(gcf, 'scatter_rf_2nd_deriv_validated.png');
% 
% % Display metrics
% fprintf('\nRandom Forest (2nd Derivatives, In-Situ Validation) Performance:\n');
% fprintf('R²: %.4f\n', r2);
% fprintf('RMSE: %.4f mg/L\n', rmse);
% fprintf('Bias: %.4f mg/L\n', mean(y_pred - tss_insitu));
% fprintf('Sample size: %d\n', length(tss_insitu));

%% train with simulated and test with insitu global data (gloria and LAB)
% Best TSS match-ups
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

<<<<<<< Updated upstream
=======
% % Load real data
% data_table = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
% tss = data_table.TSS;  % TSS values (mg/L)
% rrs_data = table2array(data_table(:, 4:end));  % Rrs spectra (rows: samples, cols: wavelengths)
% num_samples = size(rrs_data, 1);

>>>>>>> Stashed changes
% Load simulated data (training)
data_table_sim = readtable('close_tss_rrs_spectra.csv', 'VariableNamingRule', 'preserve');
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
<<<<<<< Updated upstream
% Ensure Rrs is positive
rrs_data_sim = max(rrs_data_sim, eps);

% Load in-situ data (testing)
data_table_insitu = readtable('Subset_Simulated_close_tss_rrs_spectra.csv', 'VariableNamingRule', 'preserve');
=======


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


% Load in-situ data (testing)
% data_table_insitu = readtable('Simul_Rrs_Nechad_et_al_2015_PACE_Bands_only.csv', 'VariableNamingRule', 'preserve');
% data_table_insitu = readtable('GLORIA_LAB_Rrs_02Sept2025_all_Vars_400_720TSS_PACE.csv', 'VariableNamingRule', 'preserve');
%data_table_insitu = readtable('Subset_global_GLORIA_LAB_close_tss5.csv', 'VariableNamingRule', 'preserve');
% data_table_insitu = readtable('Filtered_Subset_global_GLORIA_LAB_close_tss10.csv', 'VariableNamingRule', 'preserve');
% data_table_insitu = readtable('Subset_simulated_RRS_close_tss20.csv', 'VariableNamingRule', 'preserve');
% data_table_insitu = readtable('Filtered_Subset_global_GLORIA_LAB_close_tss10.csv', 'VariableNamingRule', 'preserve');
% data_table_insitu = readtable('GLORIA_LAB_Rrs_02Sept2025_all_Vars_400_720TSS_PACE_SAB.csv', 'VariableNamingRule', 'preserve');
data_table_insitu = readtable('Filtered_Subset_global_GLORIA_LAB_close_tss5.csv', 'VariableNamingRule', 'preserve');
>>>>>>> Stashed changes
tss_insitu = data_table_insitu.TSS;  % TSS values (mg/L)
rrs_data_insitu = table2array(data_table_insitu(:, 4:end));  % Rrs spectra
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
    try
        spec = sgolayfilt(spec, 3, 11);  % 3rd-order, 11-point window
    catch
        spec = smoothdata(spec, 'movmean', 5);  % Fallback
    end
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
best_rs = 42;
best_leaf = 4;
best_trees = 100;

% Train model on full simulated data
rng(best_rs);
rf_model = TreeBagger(best_trees, second_deriv_sim, tss_sim, ...
    'Method', 'regression', ...
    'MinLeafSize', best_leaf,...
    'OOBPredictorImportance', 'on');

% Save the model
save('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', ...
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
saveas(gcf, 'feature_importance_2nd_deriv_validated.png');

% Display feature importance summary
fprintf('\nFeature Importance Summary:\n');
fprintf('Mean importance: %.4f\n', mean(importances));
fprintf('Max importance: %.4f (at %d nm)\n', max(importances), wavelengths(imp_idx(1)));
fprintf('Min importance: %.4f (at %d nm)\n', min(importances), wavelengths(imp_idx(end)));
fprintf('Importance range: %.4f\n', range(importances));

% Model diagnostics
fprintf('\nModel Diagnostics:\n');
fprintf('Number of Observations (simulated, training): %d\n', num_samples_sim);
fprintf('Number of Observations (in-situ, testing): %d\n', num_samples_insitu);

% Scatter plot for predictions
figure(15);
hold on;
scatter(tss_insitu, y_pred, 50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
plot([min(tss_insitu), max(tss_insitu)], [min(tss_insitu), max(tss_insitu)], 'r--', 'LineWidth', 2);
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
xlabel('In-Situ TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg/L)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-2, 400]);
ylim([-2, 400]);
title('Random Forest: In-Situ vs Predicted TSS (2nd Derivatives, In-Situ Validation)', ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
p = polyfit(tss_insitu, y_pred, 1);
regress_line = polyval(p, xlim);
plot(xlim, regress_line, 'k-', 'LineWidth', 1.5, 'LineStyle', ':');
slope = p(1);
intercept = p(2);
legend({'Data Points', '1:1 Line', sprintf('Regression: y = %.2fx + %.2f', slope, intercept)}, ...
       'Location', 'southeast');
hold off;
set(gcf, 'Position', [100, 100, 800, 700]);
saveas(gcf, 'scatter_rf_2nd_deriv_validated.png');

% Display metrics
fprintf('\nRandom Forest (2nd Derivatives, In-Situ Validation) Performance:\n');
fprintf('R²: %.4f\n', r2);
fprintf('RMSE: %.4f mg/L\n', rmse);
fprintf('Bias: %.4f mg/L\n', mean(y_pred - tss_insitu));
fprintf('Sample size: %d\n', length(tss_insitu));



% %% Find samples close to 1:1 line (differences <5)
% close_idx = abs(y_pred - tss_insitu) < 5;
% num_close = sum(close_idx);
% fprintf('Found %d samples with TSS difference <5 mg/L.\n', num_close);
% 
% % Write close samples to CSV
% if num_close > 0
%     rrs_close = rrs_data_insitu(close_idx, :);
%     tss_close = tss_insitu(close_idx);
%     close_table = array2table(rrs_close, 'VariableNames', arrayfun(@(x) sprintf('Rrs_%d', x), wavelengths, 'UniformOutput', false));
%     close_table.TSS = tss_close;
%     writetable(close_table, 'Subset_Simulated_close_tss_rrs_spectra.csv');
%     fprintf('Saved close Rrs spectra and TSS to close_tss_rrs_spectra.csv\n');
% else
%     fprintf('No samples with TSS difference <5 mg/L found.\n');
% end



