%% Step 1: Set Environment, Defaults, & Paths
clc; clear; close all;

% Set graphic defaults to Arial for Q1 journal standards
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% Define Target Directory & Input File Paths
target_dir = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end
cd(target_dir);

mat_file = '/Users/masud/WorkingFile/QWIP/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_16Oct26.mat';
csv_file = fullfile(target_dir, 'St_cruise_SoRad_Matchups.csv');

% Load TreeBagger model and training wavelengths array
fprintf('Loading PACE RF model from: %s\n', mat_file);
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load SoRad Matchup CSV & Map Wavelength Columns
fprintf('Loading SoRad matchup data from: %s\n', csv_file);
matchups = readtable(csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Construct expected column names for the model's 117 wavelengths
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Check columns and fallback to nearest available Rrs column if needed
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    warning('Some model wavelengths missing in CSV; mapping to nearest available columns.');
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
    end
end

rrs_matrix = table2array(matchups(:, expected_cols));
fprintf('Extracted %d Rrs columns matching the RF model wavelengths.\n', size(rrs_matrix, 2));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

fprintf('Applying Savitzky-Golay smoothing & 2nd derivative calculation...\n');
for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21); % 3rd-order, 21-point Savitzky-Golay
end

first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Apply Original RF Model Directly (Raw Predictions)
fprintf('Predicting TSS using original PACE RF Model...\n');
tss_predicted = predict(rf_model, second_deriv_matrix);
tss_predicted = max(0, min(1000, tss_predicted)); % Clip bounds [0, 1000]

%% Step 5: Statistical Performance Analysis
valid_idx = ~isnan(tss_insitu) & ~isnan(tss_predicted);
obs = tss_insitu(valid_idx);
pred = tss_predicted(valid_idx);
N_total = length(obs);

if N_total == 0
    error('No valid matching pairs of in-situ and predicted TSS found.');
end

RMSE = sqrt(mean((pred - obs).^2));
R_mat = corrcoef(obs, pred);
Rsq = R_mat(1,2)^2;

fprintf('\nSoRad-based RF Model Evaluation Stats:\n N = %d\n R^2 = %.4f\n RMSE = %.4f mg/L\n', ...
    N_total, Rsq, RMSE);

%% Step 6: Save Output Results to CSV
output_csv_name = fullfile(target_dir, 'SoRad_RF_vs_Insitu_TSS_Matchups.csv');
results_table = table((1:N_total)', obs, pred, ...
    'VariableNames', {'Point_ID', 'InSitu_TSS', 'SoRad_RF_Predicted_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved matchup results to: %s\n', output_csv_name);

%% Step 7: Plot Q1 Journal Publication Scatter Figure (Linear Scale)
fig = figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'PaperPositionMode', 'auto');
ax = axes('Parent', fig);

hold(ax, 'on');

% Color & Line Palette
color_pts = [0.12 0.47 0.71];   % Deep Blue
color_line = [0.25 0.25 0.25];  % Dark Charcoal for 1:1 Line

% Calculate dynamic axis bounds
all_vals = [obs; pred];
min_val = max(0, floor(min(all_vals) - 4));
max_val = ceil(max(all_vals) + 4);
range_span = max_val - min_val;

% 1. Reference 1:1 Line
h_line = plot(ax, [min_val, max_val], [min_val, max_val], '--', ...
    'Color', color_line, 'LineWidth', 1.2);

% 2. Scatter Matchup Points (Circle size = 35)
h_pts = scatter(ax, obs, pred, 35, color_pts, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% 3. Add Point ID Numbers
for i = 1:N_total
    text(ax, obs(i) + (range_span * 0.018), pred(i), sprintf('%d', i), ...
        'FontSize', 8, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'middle');
end

% Formatting Axes & Grid
grid(ax, 'on');
box(ax, 'on');
set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 1.0, ...
    'XColor', 'k', 'YColor', 'k');

xlim(ax, [min_val max_val]);
ylim(ax, [min_val max_val]);

xlabel(ax, 'In Situ TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);
ylabel(ax, 'Predicted TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);

% 4. Statistical Summary Box (Bottom Right)
stats_str = { ...
    sprintf('N = %d', N_total), ...
    sprintf('R^2 = %.2f', Rsq), ...
    sprintf('RMSE = %.2f mg L^{-1}', RMSE) ...
};
text(ax, min_val + (range_span * 0.65), min_val + (range_span * 0.22), stats_str, ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1 0.9], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);

% 5. Legend (Top Left)
legend(ax, [h_pts, h_line], ...
    {'SoRad Matchups', '1:1 Line'}, ...
    'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');

hold(ax, 'off');

%% Step 8: Export Figure to Vector PDF and High-Res PNG
pdf_filename = fullfile(target_dir, 'Fig_SoRad_RF_TSS_Matchups_Oct2026.pdf');
png_filename = fullfile(target_dir, 'Fig_SoRad_RF_TSS_Matchups_Oct2026.png');

% Save Vector PDF for journal submission
exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

% Save 300 DPI PNG
exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);




%% Calibration of the model
%% Step 1: Set Environment, Defaults, & Paths
clc; clear; close all;

% Set graphic defaults to Arial for Q1 journal standards
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% Define Target Directory & Input File Paths
target_dir = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end
cd(target_dir);

mat_file = '/Users/masud/WorkingFile/QWIP/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_16Oct26.mat';
csv_file = fullfile(target_dir, 'St_cruise_SoRad_Matchups.csv');

% Load TreeBagger model and training wavelengths array
fprintf('Loading PACE RF model from: %s\n', mat_file);
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load SoRad Matchup CSV & Map Wavelength Columns
fprintf('Loading SoRad matchup data from: %s\n', csv_file);
matchups = readtable(csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Construct expected column names for the model's 117 wavelengths
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Check columns and fallback to nearest available Rrs column if needed
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    warning('Some model wavelengths missing in CSV; mapping to nearest available columns.');
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
    end
end

rrs_matrix = table2array(matchups(:, expected_cols));
fprintf('Extracted %d Rrs columns matching the RF model wavelengths.\n', size(rrs_matrix, 2));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

fprintf('Applying Savitzky-Golay smoothing & 2nd derivative calculation...\n');
for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21); % 3rd-order, 21-point Savitzky-Golay
end

first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Apply Original RF Model & Perform Regional Calibration Fit
fprintf('Predicting raw TSS using PACE RF Model...\n');
tss_raw = predict(rf_model, second_deriv_matrix);
tss_raw = tss_raw(:); % Ensure column vector format

% Valid matchup pair filter
valid_idx = ~isnan(tss_insitu) & ~isnan(tss_raw);
obs = tss_insitu(valid_idx);
pred_raw = tss_raw(valid_idx);
N_total = length(obs);

if N_total == 0
    error('No valid matching pairs of in-situ and predicted TSS found.');
end

% Perform Ordinary Least Squares Linear Calibration Fit: InSitu = slope * Raw_RF + intercept
p = polyfit(pred_raw, obs, 1);
slope = p(1);
intercept = p(2);

fprintf('\nApplied Regional Linear Calibration:\n');
fprintf('Calibrated TSS = (%.4f * Raw_SoRad_RF_TSS) + %.4f\n', slope, intercept);

% Calculate final calibrated predictions (boundary clipped at 0)
tss_predicted = max(0, (pred_raw * slope) + intercept);

%% Step 5: Statistical Performance Analysis
RMSE = sqrt(mean((tss_predicted - obs).^2));
R_mat = corrcoef(obs, tss_predicted);
Rsq = R_mat(1,2)^2;

fprintf('\nCalibrated SoRad RF Model Evaluation Stats:\n N = %d\n R^2 = %.4f\n RMSE = %.4f mg/L\n', ...
    N_total, Rsq, RMSE);

%% Step 6: Save Output Results to CSV (New Output File)
output_csv_name = fullfile(target_dir, 'SoRad_RF_Calibrated_vs_Insitu_TSS_Matchups.csv');
results_table = table((1:N_total)', obs, pred_raw, tss_predicted, ...
    'VariableNames', {'Point_ID', 'InSitu_TSS', 'Raw_SoRad_RF_TSS', 'Calibrated_SoRad_RF_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved calibrated matchup results to: %s\n', output_csv_name);

%% Step 7: Plot Q1 Journal Publication Scatter Figure (Linear Scale)
fig = figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'PaperPositionMode', 'auto');
ax = axes('Parent', fig);

hold(ax, 'on');

% Color Palette & Styling
color_pts = [0.12 0.47 0.71];   % Deep Blue
color_line = [0.25 0.25 0.25];  % Dark Charcoal for 1:1 Line

% Calculate dynamic axis bounds
all_vals = [obs; tss_predicted];
min_val = max(0, floor(min(all_vals) - 4));
max_val = ceil(max(all_vals) + 4);
range_span = max_val - min_val;

% 1. Reference 1:1 Line
h_line = plot(ax, [min_val, max_val], [min_val, max_val], '--', ...
    'Color', color_line, 'LineWidth', 1.2);

% 2. Scatter Matchup Points (Circle size = 35)
h_pts = scatter(ax, obs, tss_predicted, 35, color_pts, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% 3. Add Point ID Numbers
for i = 1:N_total
    text(ax, obs(i) + (range_span * 0.018), tss_predicted(i), sprintf('%d', i), ...
        'FontSize', 8, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'middle');
end

% Formatting Axes & Grid
grid(ax, 'on');
box(ax, 'on');
set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 1.0, ...
    'XColor', 'k', 'YColor', 'k');

xlim(ax, [min_val max_val]);
ylim(ax, [min_val max_val]);

xlabel(ax, 'In Situ TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);
ylabel(ax, 'Predicted TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);

% 4. Statistical Summary Box (Bottom Right)
stats_str = { ...
    sprintf('N = %d', N_total), ...
    sprintf('R^2 = %.2f', Rsq), ...
    sprintf('RMSE = %.2f mg L^{-1}', RMSE) ...
};
text(ax, min_val + (range_span * 0.65), min_val + (range_span * 0.22), stats_str, ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1 0.9], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);

% 5. Legend (Top Left)
legend(ax, [h_pts, h_line], ...
    {'SoRad Matchups', '1:1 Line'}, ...
    'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');

hold(ax, 'off');

%% Step 8: Export Figures to Vector PDF and High-Res PNG (New Names)
pdf_filename = fullfile(target_dir, 'Fig_SoRad_RF_Calibrated_TSS_Matchups_Oct2026.pdf');
png_filename = fullfile(target_dir, 'Fig_SoRad_RF_Calibrated_TSS_Matchups_Oct2026.png');

% Save Vector PDF for publication submission
exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

% Save 300 DPI PNG
exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);




%% Save the calibrated model.
%% Step 1: Set Environment, Defaults, & Paths
clc; clear; close all;

% Set graphic defaults to Arial for Q1 journal standards
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% Define Target Directory & Input File Paths
target_dir = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end
cd(target_dir);

mat_file = '/Users/masud/WorkingFile/QWIP/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_16Oct26.mat';
csv_file = fullfile(target_dir, 'St_cruise_SoRad_Matchups.csv');

% Load TreeBagger model and training wavelengths array
fprintf('Loading PACE RF model from: %s\n', mat_file);
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load SoRad Matchup CSV & Map Wavelength Columns
fprintf('Loading SoRad matchup data from: %s\n', csv_file);
matchups = readtable(csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Construct expected column names for the model's 117 wavelengths
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Check columns and fallback to nearest available Rrs column if needed
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    warning('Some model wavelengths missing in CSV; mapping to nearest available columns.');
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
    end
end

rrs_matrix = table2array(matchups(:, expected_cols));
fprintf('Extracted %d Rrs columns matching the RF model wavelengths.\n', size(rrs_matrix, 2));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

fprintf('Applying Savitzky-Golay smoothing & 2nd derivative calculation...\n');
for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21); % 3rd-order, 21-point Savitzky-Golay
end

first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Apply Original RF Model & Perform Regional Calibration Fit
fprintf('Predicting raw TSS using PACE RF Model...\n');
tss_raw = predict(rf_model, second_deriv_matrix);
tss_raw = tss_raw(:); % Ensure column vector format

% Valid matchup pair filter
valid_idx = ~isnan(tss_insitu) & ~isnan(tss_raw);
obs = tss_insitu(valid_idx);
pred_raw = tss_raw(valid_idx);
N_total = length(obs);

if N_total == 0
    error('No valid matching pairs of in-situ and predicted TSS found.');
end

% Perform Ordinary Least Squares Linear Calibration Fit: InSitu = slope * Raw_RF + intercept
p = polyfit(pred_raw, obs, 1);
slope = p(1);
intercept = p(2);

fprintf('\nApplied Regional Linear Calibration:\n');
fprintf('Calibrated TSS = (%.4f * Raw_SoRad_RF_TSS) + %.4f\n', slope, intercept);

% Calculate final calibrated predictions (boundary clipped at 0)
tss_predicted = max(0, (pred_raw * slope) + intercept);

%% Step 5: Save Calibrated Model Object (.mat File)
calibrated_mat_file = fullfile(target_dir, 'PACE_rf_model_2nd_deriv_SoRad_Calibrated_Oct2026.mat');

% Build clean model structure container
calibrated_rf_model = struct();
calibrated_rf_model.rf_model = rf_model;
calibrated_rf_model.wavelengths = wavelengths;
calibrated_rf_model.calibration_slope = slope;
calibrated_rf_model.calibration_intercept = intercept;
calibrated_rf_model.equation = sprintf('TSS_calibrated = max(0, %.6f * TSS_raw + %.6f)', slope, intercept);
calibrated_rf_model.calibration_date = datestr(now);

% Save to disk
save(calibrated_mat_file, 'calibrated_rf_model', 'rf_model', 'wavelengths', 'slope', 'intercept');
fprintf('\nSuccessfully saved calibrated model structure to:\n%s\n', calibrated_mat_file);

%% Step 6: Statistical Performance Analysis
RMSE = sqrt(mean((tss_predicted - obs).^2));
R_mat = corrcoef(obs, tss_predicted);
Rsq = R_mat(1,2)^2;

fprintf('\nCalibrated SoRad RF Model Evaluation Stats:\n N = %d\n R^2 = %.4f\n RMSE = %.4f mg/L\n', ...
    N_total, Rsq, RMSE);

%% Step 7: Save Output Results to CSV
output_csv_name = fullfile(target_dir, 'SoRad_RF_Calibrated_vs_Insitu_TSS_Matchups.csv');
results_table = table((1:N_total)', obs, pred_raw, tss_predicted, ...
    'VariableNames', {'Point_ID', 'InSitu_TSS', 'Raw_SoRad_RF_TSS', 'Calibrated_SoRad_RF_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved calibrated matchup results to: %s\n', output_csv_name);

%% Step 8: Plot Q1 Journal Publication Scatter Figure (Linear Scale)
fig = figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'PaperPositionMode', 'auto');
ax = axes('Parent', fig);

hold(ax, 'on');

% Color Palette & Styling
color_pts = [0.12 0.47 0.71];   % Deep Blue
color_line = [0.25 0.25 0.25];  % Dark Charcoal for 1:1 Line

% Calculate dynamic axis bounds
all_vals = [obs; tss_predicted];
min_val = max(0, floor(min(all_vals) - 4));
max_val = ceil(max(all_vals) + 4);
range_span = max_val - min_val;

% 1. Reference 1:1 Line
h_line = plot(ax, [min_val, max_val], [min_val, max_val], '--', ...
    'Color', color_line, 'LineWidth', 1.2);

% 2. Scatter Matchup Points (Circle size = 35)
h_pts = scatter(ax, obs, tss_predicted, 35, color_pts, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% 3. Add Point ID Numbers
for i = 1:N_total
    text(ax, obs(i) + (range_span * 0.018), tss_predicted(i), sprintf('%d', i), ...
        'FontSize', 8, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'middle');
end

% Formatting Axes & Grid
grid(ax, 'on');
box(ax, 'on');
set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 1.0, ...
    'XColor', 'k', 'YColor', 'k');

xlim(ax, [min_val max_val]);
ylim(ax, [min_val max_val]);

xlabel(ax, 'In Situ TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);
ylabel(ax, 'Predicted TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);

% 4. Statistical Summary Box (Bottom Right)
stats_str = { ...
    sprintf('N = %d', N_total), ...
    sprintf('R^2 = %.2f', Rsq), ...
    sprintf('RMSE = %.2f mg L^{-1}', RMSE) ...
};
text(ax, min_val + (range_span * 0.65), min_val + (range_span * 0.22), stats_str, ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1 0.9], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);

% 5. Legend (Top Left)
legend(ax, [h_pts, h_line], ...
    {'SoRad Matchups', '1:1 Line'}, ...
    'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');

hold(ax, 'off');

%% Step 9: Export Figures to Vector PDF and High-Res PNG
pdf_filename = fullfile(target_dir, 'Fig_SoRad_RF_Calibrated_TSS_Matchups_Oct2026.pdf');
png_filename = fullfile(target_dir, 'Fig_SoRad_RF_Calibrated_TSS_Matchups_Oct2026.png');

% Save Vector PDF for publication submission
exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

% Save 300 DPI PNG
exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);


%% Apply calibrated model to PACE
%% Step 1: Set Environment, Graphic Defaults, & File Paths
clc; clear; close all;

% Set publication graphic defaults (Arial font for Q1 journals)
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% File Paths
sorad_dir = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';
calibrated_mat_file = fullfile(sorad_dir, 'PACE_rf_model_2nd_deriv_SoRad_Calibrated_Oct2026.mat');
pace_csv_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/CleanedRrs/St_cruise_PACE_26_09_2026.csv';

% Load Calibrated RF Model Object & Parameters
fprintf('Loading Calibrated PACE RF Model from: %s\n', calibrated_mat_file);
load(calibrated_mat_file, 'calibrated_rf_model');

% Extract components from calibrated structure
rf_model  = calibrated_rf_model.rf_model;
model_wl  = calibrated_rf_model.wavelengths;
slope     = calibrated_rf_model.calibration_slope;
intercept = calibrated_rf_model.calibration_intercept;

fprintf('Calibration Equation: TSS_calibrated = max(0, %.4f * TSS_raw + %.4f)\n', slope, intercept);

%% Step 2: Load PACE Matchup CSV & Extract Wavelength Columns
fprintf('Loading PACE matchup data from: %s\n', pace_csv_file);
matchups = readtable(pace_csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Construct expected column names for the model's 117 wavelengths
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Check columns and fallback to nearest available Rrs column if needed
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    warning('Some model wavelengths missing in CSV; mapping to nearest available columns.');
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
    end
end

rrs_matrix = table2array(matchups(:, expected_cols));
fprintf('Extracted %d Rrs columns matching the RF model wavelengths.\n', size(rrs_matrix, 2));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

fprintf('Applying Savitzky-Golay smoothing & 2nd derivative calculation...\n');
for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21); % 3rd-order, 21-point Savitzky-Golay
end

first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Apply Calibrated Model (Raw Prediction + Slope & Intercept Correction)
fprintf('Predicting TSS using Random Forest model...\n');
tss_raw = predict(rf_model, second_deriv_matrix);
tss_raw = tss_raw(:); % Ensure column vector format

% Apply calibration transformation
tss_predicted = max(0, (tss_raw * slope) + intercept);

%% Step 5: Statistical Performance Analysis
valid_idx = ~isnan(tss_insitu) & ~isnan(tss_predicted);
obs = tss_insitu(valid_idx);
pred = tss_predicted(valid_idx);
pred_raw = tss_raw(valid_idx);
N_total = length(obs);

if N_total == 0
    error('No valid matching pairs of in-situ and predicted TSS found.');
end

RMSE = sqrt(mean((pred - obs).^2));
R_mat = corrcoef(obs, pred);
Rsq = R_mat(1,2)^2;

fprintf('\nCalibrated PACE Model Evaluation Stats:\n N = %d\n R^2 = %.4f\n RMSE = %.4f mg/L\n', ...
    N_total, Rsq, RMSE);

%% Step 6: Save Output Results to CSV
output_csv_name = 'Calibrated_RF_vs_Insitu_PACE_TSS_Matchups.csv';
results_table = table((1:N_total)', obs, pred_raw, pred, ...
    'VariableNames', {'Point_ID', 'InSitu_TSS', 'Raw_PACE_RF_TSS', 'Calibrated_PACE_RF_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved calibrated PACE matchup results to: %s\n', output_csv_name);

%% Step 7: Plot Q1 Journal Publication Scatter Figure (Linear Scale)
fig = figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'PaperPositionMode', 'auto');
ax = axes('Parent', fig);

hold(ax, 'on');

% Colors & Styling
color_pts = [0.12 0.47 0.71];   % Deep Blue
color_line = [0.25 0.25 0.25];  % Dark Charcoal for 1:1 Line

% Calculate dynamic axis bounds
all_vals = [obs; pred];
min_val = max(0, floor(min(all_vals) - 4));
max_val = ceil(max(all_vals) + 4);
range_span = max_val - min_val;

% 1. Reference 1:1 Line
h_line = plot(ax, [min_val, max_val], [min_val, max_val], '--', ...
    'Color', color_line, 'LineWidth', 1.2);

% 2. Scatter Matchup Points (Circle size = 35)
h_pts = scatter(ax, obs, pred, 35, color_pts, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% 3. Add Point ID Numbers
for i = 1:N_total
    text(ax, obs(i) + (range_span * 0.018), pred(i), sprintf('%d', i), ...
        'FontSize', 8, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'middle');
end

% Formatting Axes & Grid
grid(ax, 'on');
box(ax, 'on');
set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 1.0, ...
    'XColor', 'k', 'YColor', 'k');

xlim(ax, [min_val max_val]);
ylim(ax, [min_val max_val]);

xlabel(ax, 'In Situ TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);
ylabel(ax, 'Predicted TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);

% 4. Statistical Summary Box (Bottom Right)
stats_str = { ...
    sprintf('N = %d', N_total), ...
    sprintf('R^2 = %.2f', Rsq), ...
    sprintf('RMSE = %.2f mg L^{-1}', RMSE) ...
};
text(ax, min_val + (range_span * 0.65), min_val + (range_span * 0.22), stats_str, ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1 0.9], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);

% 5. Legend (Top Left)
legend(ax, [h_pts, h_line], ...
    {'PACE Matchups', '1:1 Line'}, ...
    'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');

hold(ax, 'off');

%% Step 8: Export Figures to Vector PDF and High-Res PNG
pdf_filename = 'Fig_Calibrated_RF_PACE_TSS_Matchups.pdf';
png_filename = 'Fig_Calibrated_RF_PACE_TSS_Matchups.png';

% Save Vector PDF for journal submission
exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

% Save 300 DPI PNG
exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);


% %% Sorad vs PACE spectra:
% %% Step 1: Initialize Workspace & Set Publication Typography Defaults
% clc; clear; close all;
% 
% % Set publication graphic defaults (Arial font for Q1 journals)
% set(0, 'DefaultAxesFontName', 'Arial');
% set(0, 'DefaultTextFontName', 'Arial');
% set(0, 'DefaultAxesFontSize', 10);
% set(0, 'DefaultTextFontSize', 10);
% 
% % File Paths
% pace_file  = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/CleanedRrs/St_cruise_PACE_26_09_2026.csv';
% sorad_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad/St_cruise_SoRad_Matchups.csv';
% 
% %% Step 2: Read CSV Files and Parse Wavelength Headers
% fprintf('Loading PACE data: %s\n', pace_file);
% pace_table = readtable(pace_file);
% 
% fprintf('Loading SoRad data: %s\n', sorad_file);
% sorad_table = readtable(sorad_file);
% 
% % Extract PACE Rrs column names and numeric wavelengths
% pace_vars = pace_table.Properties.VariableNames;
% pace_rrs_cols = pace_vars(startsWith(pace_vars, 'Rrs_', 'IgnoreCase', true));
% pace_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), pace_rrs_cols);
% 
% % Extract SoRad Rrs column names and numeric wavelengths
% sorad_vars = sorad_table.Properties.VariableNames;
% sorad_rrs_cols = sorad_vars(startsWith(sorad_vars, 'Rrs_', 'IgnoreCase', true));
% sorad_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), sorad_rrs_cols);
% 
% %% Step 3: Isolate First 7 Stations Spectra
% N_plot = min([7, height(pace_table), height(sorad_table)]);
% 
% pace_spectra  = table2array(pace_table(1:N_plot, pace_rrs_cols));
% sorad_spectra = table2array(sorad_table(1:N_plot, sorad_rrs_cols));
% 
% %% Step 4: Create Publication Figure Overlaid Plot
% fig = figure('Units', 'centimeters', 'Position', [3, 3, 14, 10], 'PaperPositionMode', 'auto');
% ax = axes('Parent', fig);
% hold(ax, 'on');
% 
% % Define 7 distinct color palette (Q1 colorblind friendly)
% color_palette = [
%     0.12, 0.47, 0.71;  % St 1: Blue
%     0.89, 0.10, 0.11;  % St 2: Red
%     0.17, 0.63, 0.17;  % St 3: Green
%     1.00, 0.50, 0.00;  % St 4: Orange
%     0.58, 0.40, 0.74;  % St 5: Purple
%     0.65, 0.34, 0.16;  % St 6: Brown
%     0.88, 0.47, 0.76   % St 7: Pink
% ];
% 
% h_pace = zeros(1, N_plot);
% h_sorad = zeros(1, N_plot);
% 
% % Plot spectra for each station
% for i = 1:N_plot
%     c = color_palette(i, :);
% 
%     % PACE spectrum (Solid line, width = 1.4)
%     h_pace(i) = plot(ax, pace_wls, pace_spectra(i, :), '-', ...
%         'Color', c, 'LineWidth', 1.4);
% 
%     % SoRad spectrum (Dashed line, width = 1.2)
%     h_sorad(i) = plot(ax, sorad_wls, sorad_spectra(i, :), '--', ...
%         'Color', c, 'LineWidth', 1.2);
% end
% 
% % Formatting Axes, Grid, & Labels
% grid(ax, 'on');
% box(ax, 'on');
% set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 1.0, ...
%     'XColor', 'k', 'YColor', 'k');
% 
% xlim(ax, [340 730]);
% xlabel(ax, 'Wavelength \lambda (nm)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel(ax, 'Remote Sensing Reflectance R_{rs}(\lambda) (sr^{-1})', 'FontWeight', 'bold', 'FontSize', 10);
% 
% % Build Clean Custom Legend
% legend_handles = [];
% legend_labels = {};
% 
% for i = 1:N_plot
%     legend_handles = [legend_handles, h_pace(i), h_sorad(i)]; %#ok<AGROW>
%     legend_labels = [legend_labels, {sprintf('St %d (PACE)', i), sprintf('St %d (SoRad)', i)}]; %#ok<AGROW>
% end
% 
% legend(ax, legend_handles, legend_labels, ...
%     'NumColumns', 1, 'Location', 'northeastoutside', ...
%     'FontSize', 8, 'Box', 'on');
% 
% hold(ax, 'off');
% 
% %% Step 5: Export Vector PDF & High-Res PNG
% pdf_filename = 'Fig_PACE_vs_SoRad_Spectra_Q1.pdf';
% png_filename = 'Fig_PACE_vs_SoRad_Spectra_Q1.png';
% 
% exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
% fprintf('Exported Vector PDF: %s\n', pdf_filename);
% 
% exportgraphics(fig, png_filename, 'Resolution', 300);
% fprintf('Exported 300 DPI PNG: %s\n', png_filename);

% %% Station wise
% %% Step 1: Initialize Workspace & Set Publication Typography Defaults
% clc; clear; close all;
% 
% % Set publication graphic defaults (Arial font for Q1 journals)
% set(0, 'DefaultAxesFontName', 'Arial');
% set(0, 'DefaultTextFontName', 'Arial');
% set(0, 'DefaultAxesFontSize', 9);
% set(0, 'DefaultTextFontSize', 9);
% 
% % File Paths
% pace_file  = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/CleanedRrs/St_cruise_PACE_26_09_2026.csv';
% sorad_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad/St_cruise_SoRad_Matchups_Filtered_10min.csv';
% 
% %% Step 2: Read CSV Files and Parse Wavelength Headers
% fprintf('Loading PACE data: %s\n', pace_file);
% pace_table = readtable(pace_file);
% 
% fprintf('Loading SoRad data: %s\n', sorad_file);
% sorad_table = readtable(sorad_file);
% 
% % Extract PACE Rrs column names and numeric wavelengths
% pace_vars = pace_table.Properties.VariableNames;
% pace_rrs_cols = pace_vars(startsWith(pace_vars, 'Rrs_', 'IgnoreCase', true));
% pace_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), pace_rrs_cols);
% 
% % Extract SoRad Rrs column names and numeric wavelengths
% sorad_vars = sorad_table.Properties.VariableNames;
% sorad_rrs_cols = sorad_vars(startsWith(sorad_vars, 'Rrs_', 'IgnoreCase', true));
% sorad_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), sorad_rrs_cols);
% 
% %% Step 3: Filter Wavelengths (400 to 719 nm) & Isolate First 7 Stations
% % Wavelength masks for 400 - 719 nm
% pace_wl_mask  = (pace_wls >= 400) & (pace_wls <= 719);
% sorad_wl_mask = (sorad_wls >= 400) & (sorad_wls <= 719);
% 
% pace_wls_sub  = pace_wls(pace_wl_mask);
% sorad_wls_sub = sorad_wls(sorad_wl_mask);
% 
% N_plot = min([7, height(pace_table), height(sorad_table)]);
% 
% pace_spectra  = table2array(pace_table(1:N_plot, pace_rrs_cols));
% sorad_spectra = table2array(sorad_table(1:N_plot, sorad_rrs_cols));
% 
% pace_spectra_sub  = pace_spectra(:, pace_wl_mask);
% sorad_spectra_sub = sorad_spectra(:, sorad_wl_mask);
% 
% %% Step 4: Create Multi-Panel Subplot Figure (4x2 Grid)
% fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 22], 'PaperPositionMode', 'auto');
% 
% % Define 7 distinct Q1 colorblind-friendly colors
% color_palette = [
%     0.12, 0.47, 0.71;  % St 1: Blue
%     0.89, 0.10, 0.11;  % St 2: Red
%     0.17, 0.63, 0.17;  % St 3: Green
%     0.80, 0.40, 0.00;  % St 4: Dark Orange
%     0.58, 0.40, 0.74;  % St 5: Purple
%     0.65, 0.34, 0.16;  % St 6: Brown
%     0.88, 0.47, 0.76   % St 7: Pink
% ];
% 
% for i = 1:N_plot
%     ax = subplot(4, 2, i);
%     hold(ax, 'on');
% 
%     c = color_palette(i, :);
% 
%     % Plot PACE spectrum (Solid Line)
%     h_p = plot(ax, pace_wls_sub, pace_spectra_sub(i, :), '-', ...
%         'Color', c, 'LineWidth', 1.4);
% 
%     % Plot SoRad spectrum (Dashed Line)
%     h_s = plot(ax, sorad_wls_sub, sorad_spectra_sub(i, :), '--', ...
%         'Color', c * 0.7, 'LineWidth', 1.2);
% 
%     % Formatting Individual Axes
%     grid(ax, 'on');
%     box(ax, 'on');
%     set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 0.8, ...
%         'XColor', 'k', 'YColor', 'k');
% 
%     xlim(ax, [400 719]);
%     title(ax, sprintf('Station %d', i), 'FontWeight', 'bold', 'FontSize', 9.5);
% 
%     % Add Labels conditionally for clean layout
%     if i >= 6
%         xlabel(ax, 'Wavelength \lambda (nm)', 'FontWeight', 'bold', 'FontSize', 8.5);
%     end
%     if mod(i, 2) == 1
%         ylabel(ax, 'R_{rs} (sr^{-1})', 'FontWeight', 'bold', 'FontSize', 8.5);
%     end
% 
%     % Legend per Subplot
%     legend(ax, [h_p, h_s], {'PACE', 'SoRad'}, ...
%         'Location', 'northeast', 'FontSize', 7.5, 'Box', 'on');
% 
%     hold(ax, 'off');
% end
% 
% % Hide unused 8th subplot panel
% ax_unused = subplot(4, 2, 8);
% axis(ax_unused, 'off');
% 
% %% Step 5: Export Vector PDF & High-Res PNG
% pdf_filename = 'Fig_PACE_vs_SoRad_Spectra_Subplots_Q1_10Mins.pdf';
% png_filename = 'Fig_PACE_vs_SoRad_Spectra_Subplots_Q1_10Mins.png';
% 
% exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
% fprintf('Exported Vector PDF: %s\n', pdf_filename);
% 
% exportgraphics(fig, png_filename, 'Resolution', 300);
% fprintf('Exported 300 DPI PNG: %s\n', png_filename);

%% Time stamp
%% Step 1: Initialize Workspace & Set Publication Typography Defaults
clc; clear; close all;

% Set publication graphic defaults (Arial font for Q1 journals)
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 9);
set(0, 'DefaultTextFontSize', 9);

% Define Working Directory & File Paths
sorad_dir  = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';
pace_file  = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/CleanedRrs/St_cruise_PACE_26_09_2026.csv';
sorad_file = fullfile(sorad_dir, 'St_cruise_SoRad_Matchups_Filtered_10min.csv');

%% Step 2: Read CSV Files and Parse Wavelength Headers
fprintf('Loading PACE data: %s\n', pace_file);
pace_table = readtable(pace_file);

fprintf('Loading SoRad data: %s\n', sorad_file);
sorad_table = readtable(sorad_file);

% Extract PACE Rrs column names and numeric wavelengths
pace_vars = pace_table.Properties.VariableNames;
pace_rrs_cols = pace_vars(startsWith(pace_vars, 'Rrs_', 'IgnoreCase', true));
pace_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), pace_rrs_cols);

% Extract SoRad Rrs column names and numeric wavelengths
sorad_vars = sorad_table.Properties.VariableNames;
sorad_rrs_cols = sorad_vars(startsWith(sorad_vars, 'Rrs_', 'IgnoreCase', true));
sorad_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), sorad_rrs_cols);

%% Step 3: Filter Wavelengths (400 to 719 nm) & Isolate First 7 Stations
pace_wl_mask  = (pace_wls >= 400) & (pace_wls <= 719);
sorad_wl_mask = (sorad_wls >= 400) & (sorad_wls <= 719);

pace_wls_sub  = pace_wls(pace_wl_mask);
sorad_wls_sub = sorad_wls(sorad_wl_mask);

N_plot = min([7, height(pace_table), height(sorad_table)]);

pace_spectra  = table2array(pace_table(1:N_plot, pace_rrs_cols));
sorad_spectra = table2array(sorad_table(1:N_plot, sorad_rrs_cols));

pace_spectra_sub  = pace_spectra(:, pace_wl_mask);
sorad_spectra_sub = sorad_spectra(:, sorad_wl_mask);

%% Step 4: Parse Timestamps for PACE and SoRad
% SoRad Timestamps (EDT Local Time)
if ismember('SoRad_TimeEST', sorad_table.Properties.VariableNames)
    sorad_time_str_cell = cellstr(string(sorad_table.SoRad_TimeEST(1:N_plot)));
else
    sorad_time_str_cell = cellstr(string(sorad_table.TimeEST(1:N_plot)));
end

% Station Dates
if ismember('Date', sorad_table.Properties.VariableNames)
    date_str_cell = cellstr(string(sorad_table.Date(1:N_plot)));
else
    date_str_cell = {'9/26/26', '9/26/26', '9/26/26', '9/26/26', '9/27/26', '9/27/26', '9/27/26'};
end

% PACE Fixed UTC Times by Cruise Date
% Sept 26, 2026 -> 17:19:24 UTC | Sept 27, 2026 -> 17:53:57 UTC
pace_time_str_cell = cell(N_plot, 1);
for i = 1:N_plot
    d_str = date_str_cell{i};
    if contains(d_str, '9/26') || contains(d_str, '26_09') || i <= 4
        pace_time_str_cell{i} = '17:19:24 UTC';
    else
        pace_time_str_cell{i} = '17:53:57 UTC';
    end
end

%% Step 5: Create Multi-Panel Subplot Figure (4x2 Grid)
fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 22], 'PaperPositionMode', 'auto');

% Define 7 distinct Q1 colorblind-friendly colors
color_palette = [
    0.12, 0.47, 0.71;  % St 1: Blue
    0.89, 0.10, 0.11;  % St 2: Red
    0.17, 0.63, 0.17;  % St 3: Green
    0.80, 0.40, 0.00;  % St 4: Dark Orange
    0.58, 0.40, 0.74;  % St 5: Purple
    0.65, 0.34, 0.16;  % St 6: Brown
    0.88, 0.47, 0.76   % St 7: Pink
];

for i = 1:N_plot
    ax = subplot(4, 2, i);
    hold(ax, 'on');
    
    c = color_palette(i, :);
    
    % Plot PACE spectrum (Solid Line)
    h_p = plot(ax, pace_wls_sub, pace_spectra_sub(i, :), '-', ...
        'Color', c, 'LineWidth', 1.4);
    
    % Plot SoRad spectrum (Dashed Line)
    h_s = plot(ax, sorad_wls_sub, sorad_spectra_sub(i, :), '--', ...
        'Color', c * 0.7, 'LineWidth', 1.2);
    
    % Formatting Individual Axes
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 0.8, ...
        'XColor', 'k', 'YColor', 'k');
    
    xlim(ax, [400 719]);
    
    % Title with Station ID and Date
    title_str = sprintf('Station %d (%s)', i, date_str_cell{i});
    title(ax, title_str, 'FontWeight', 'bold', 'FontSize', 9);
    
    % Conditional Axis Labels
    if i >= 6
        xlabel(ax, 'Wavelength \lambda (nm)', 'FontWeight', 'bold', 'FontSize', 8.5);
    end
    if mod(i, 2) == 1
        ylabel(ax, 'R_{rs} (\lambda) (sr^{-1})', 'FontWeight', 'bold', 'FontSize', 8.5);
    end
    
    % Legend displaying exact acquisition timestamps
    pace_legend_str  = sprintf('PACE (%s)', pace_time_str_cell{i});
    sorad_legend_str = sprintf('SoRad (%s EDT)', sorad_time_str_cell{i});
    
    legend(ax, [h_p, h_s], {pace_legend_str, sorad_legend_str}, ...
        'Location', 'northeast', 'FontSize', 7.0, 'Box', 'on');
    
    hold(ax, 'off');
end

% Hide unused 8th subplot panel
ax_unused = subplot(4, 2, 8);
axis(ax_unused, 'off');

%% Step 6: Export Vector PDF & High-Res PNG
pdf_filename = fullfile(sorad_dir, 'Fig_PACE_vs_SoRad_Spectra_Subplots_Q1_10Mins_time.pdf');
png_filename = fullfile(sorad_dir, 'Fig_PACE_vs_SoRad_Spectra_Subplots_Q1_10Mins_Time.png');

exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);

