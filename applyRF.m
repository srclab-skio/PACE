% MATLAB script to apply the best Random Forest Model to the NC file for TSS prediction
% Author: Grok (xAI) - Expert in Ocean Optics and Machine Learning
% Date: August 15, 2025
% Assumptions: 'subset_PACE.nc' and 'rf_model_best.mat' in working directory
% Focus: Apply saved RF model on second derivatives from NC Rrs bands
% Updates: Fixed scaling issue; added extensive debug prints; double precision; m_map with high-res coastline; read predicted_tss.nc
% Requirements: m_map toolbox, Statistics and Machine Learning Toolbox, Signal Processing Toolbox (optional)
% MATLAB Version: R2022a

clear; close all; clc;

% Step 1: Load the saved model
load('PACE_rf_model_best_Simulated.mat', 'rf_model_best', 'wavelengths');


% Model wavelengths (403 to 718 nm, 117 bands)
model_wl = wavelengths;

% Step 2: Inspect NetCDF file to get Rrs variables and their wavelengths
info = ncinfo('subset_PACE.nc');
var_names = {info.Variables.Name};

% Find Rrs variables (starting with 'Rrs_', excluding 'Rrs_unc_')
rrs_vars = var_names(startsWith(var_names, 'Rrs_') & ~startsWith(var_names, 'Rrs_unc_'));
num_rrs = length(rrs_vars);

% Extract actual wavelengths from attributes
nc_wl = zeros(1, num_rrs);
for i = 1:num_rrs
    nc_wl(i) = ncreadatt('subset_PACE.nc', rrs_vars{i}, 'radiation_wavelength');
end

% Sort by wavelength
[sorted_nc_wl, sort_idx] = sort(nc_wl);
sorted_rrs_vars = rrs_vars(sort_idx);

% Step 3: Map model wavelengths to nearest NC wavelengths/variables
mapped_vars = cell(1, length(model_wl));
mapped_nc_wl = zeros(1, length(model_wl));
for i = 1:length(model_wl)
    [~, idx] = min(abs(sorted_nc_wl - model_wl(i)));
    mapped_vars{i} = sorted_rrs_vars{idx};
    mapped_nc_wl(i) = sorted_nc_wl(idx);
%     fprintf('Model wl %d mapped to NC wl %.3f (%s)\n', model_wl(i), mapped_nc_wl(i), mapped_vars{i});
end

% Step 4: Read dimensions and lat/lon
x_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'x')).Length;  % 273
y_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'y')).Length;  % 426
lon = ncread('subset_PACE.nc', 'longitude');
lat = ncread('subset_PACE.nc', 'latitude');

% Compute bounding box with padding (bigger range)
lon_min = min(lon(:)) - 0.5;  % Pad by 0.5 degrees
lon_max = max(lon(:)) + 0.5;
lat_min = min(lat(:)) - 0.5;
lat_max = max(lat(:)) + 0.5;

% Step 5: Read raw Rrs data for mapped variables into 3D array (x, y, wl)
rrs_3d = zeros(x_dim, y_dim, length(mapped_vars), 'double');  % Double for precision
for i = 1:length(mapped_vars)
    var = mapped_vars{i};
    rrs_raw = ncread('subset_PACE.nc', var);
    % Apply fill value, no scaling
    fill = ncreadatt('subset_PACE.nc', var, '_FillValue');
    rrs_raw(rrs_raw == fill) = NaN;
    rrs_3d(:,:,i) = double(rrs_raw);  % Use raw values directly
    % Debug: Print stats for each band
%     fprintf('Band %s (wl %.3f): raw min %.6f, max %.6f, std %.6f\n', ...
%         var, mapped_nc_wl(i), min(rrs_raw(:),'omitnan'), max(rrs_raw(:),'omitnan'), ...
%         std(rrs_raw(:),'omitnan'));
end
%% Scaled one:to test

% % Step 5: Read raw Rrs data for mapped variables into 3D array (x, y, wl)
% rrs_3d = zeros(x_dim, y_dim, length(mapped_vars), 'double');  % Double for precision
% for i = 1:length(mapped_vars)
%     var = mapped_vars{i};
%     rrs_raw = ncread('subset_PACE.nc', var);
%     % Apply fill value
%     fill = ncreadatt('subset_PACE.nc', var, '_FillValue');
%     rrs_raw(rrs_raw == fill) = NaN;
%     % Apply scaling and offset
%     scale = ncreadatt('subset_PACE.nc', var, 'scale_factor');
%     offset = ncreadatt('subset_PACE.nc', var, 'add_offset');
%     rrs_scaled = double(rrs_raw) * scale + offset;
%     rrs_3d(:,:,i) = rrs_scaled;
% end

%%
% Debug: Print Rrs for sample pixels (first 3 good pixels)
nan_mask = any(isnan(rrs_3d), 3);  % NaN in any band
good_pixels = find(~nan_mask);
fprintf('Sample Rrs for first 3 good pixels (wl 402.654, 555.044, 716.817):\n');
for g = 1:min(3, length(good_pixels))
    [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
%     fprintf('Pixel (%d,%d): Rrs_403=%.6f, Rrs_555=%.6f, Rrs_717=%.6f\n', ...
%         ix, iy, rrs_3d(ix,iy,1), rrs_3d(ix,iy,find(mapped_nc_wl >= 555,1)), rrs_3d(ix,iy,end));
end

% Check for l2_flags to improve masking % don't apply this for my data
% try
%     l2_flags = ncread('subset_PACE.nc', 'l2_flags');
%     invalid_mask = bitand(l2_flags, 1+2+4) > 0;  % Example: land (1), cloud (2), fail (4)
%     nan_mask = nan_mask | invalid_mask;
%     disp('Applied l2_flags masking for land/clouds');
% catch
%     disp('l2_flags not found; using only Rrs NaN mask');
% end

% Step 6: Optimization - Identify good (non-masked) pixels
good_pixels = find(~nan_mask);
num_good = length(good_pixels);
if num_good == 0
    error('No valid pixels found; all masked as invalid.');
end
fprintf('Computing for %d valid pixels (%.1f%% of total).\n', num_good, 100*num_good/(x_dim*y_dim));

% Extract Rrs for good pixels only
rrs_2d_good = zeros(num_good, length(model_wl), 'double');
for g = 1:num_good
    [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
    rrs_2d_good(g, :) = squeeze(rrs_3d(ix, iy, :))';
end

% Debug: Check variance in Rrs and plot multiple spectra
rrs_std_pixels = std(rrs_2d_good, 1, 1);
rrs_std_wl = std(rrs_2d_good, 1, 2);
fprintf('Rrs std across pixels: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_pixels), max(rrs_std_pixels), mean(rrs_std_pixels));
fprintf('Rrs std across wl: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_wl), max(rrs_std_wl), mean(rrs_std_wl));
% Plot spectra for first 3 good pixels
figure(3);
hold on;
colors = {'b-', 'r-', 'g-'};
for g = 1:min(3, num_good)
    plot(model_wl, rrs_2d_good(g,:), colors{g}, 'LineWidth', 2, 'DisplayName', sprintf('Pixel %d', g));
end
xlabel('Wavelength (nm)');
ylabel('Rrs (raw units)');
title('Sample Rrs Spectra (First 3 Good Pixels)');
legend('show');
grid on;
saveas(gcf, 'Simulated_sample_rrs_spectra.png');

% Step 7: Compute second derivatives for good pixels
second_deriv_good = zeros(num_good, length(model_wl), 'double');
for g = 1:num_good
    spec = rrs_2d_good(g, :);
    % Ensure positive values
    spec = max(spec, eps);
    % Smooth
    try
        spec = sgolayfilt(spec, 3, 11);
    catch
        spec = smoothdata(spec, 'movmean', 5);
    end
    % Derivatives
    first_deriv = gradient(spec, model_wl);
    second_deriv_good(g, :) = gradient(first_deriv, model_wl);
end

% Debug: Check variance in second derivatives
deriv_std_pixels = std(second_deriv_good, 1, 1);
deriv_std_wl = std(second_deriv_good, 1, 2);
fprintf('Second deriv std across pixels: min %.6f, max %.6f, mean %.6f\n', min(deriv_std_pixels), max(deriv_std_pixels), mean(deriv_std_pixels));
fprintf('Second deriv std across wl: min %.6f, max %.6f, mean %.6f\n', min(deriv_std_wl), max(deriv_std_wl), mean(deriv_std_wl));
fprintf('Sample second deriv for first pixel: min %.6f, max %.6f\n', min(second_deriv_good(1,:)), max(second_deriv_good(1,:)));

% Step 8: Predict TSS for good pixels
tss_pred_good = predict(rf_model_best, second_deriv_good);

% Debug: Check predictions before clip
fprintf('Predicted TSS before clip min: %.4f, max: %.4f, mean: %.4f, std: %.4f\n', min(tss_pred_good), max(tss_pred_good), mean(tss_pred_good), std(tss_pred_good));

% Create full TSS map with NaNs
tss_map = NaN(x_dim, y_dim);
tss_map(good_pixels) = tss_pred_good;

% Clip to realistic range (0-100 mg/L) only for non-masked
tss_map(~isnan(tss_map)) = max(0, min(100, tss_map(~isnan(tss_map))));

% Step 9: Spatial map plot using m_map with Lambert projection and high-res coastline
m_proj('lambert', 'long', [lon_min lon_max], 'lat', [lat_min lat_max]);
figure(1);
m_pcolor(lon, lat, tss_map);
shading flat;
hold on;
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');  % High-resolution coastline and internal waters
hold off;
h = colorbar;
h.Label.String = 'TSS (mg/L)';
h.Label.FontSize = 12;
caxis([0 50]);  % Adjusted range for coastal TSS
m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
title('Predicted TSS (mg/L) from RF Model', 'FontSize', 14);
colormap(jet);
saveas(gcf, 'Simulated_tss_map.png');

% Step 10: Save predicted TSS to new NetCDF file
new_nc = 'Simulated_predicted_tss.nc';
copyfile('subset_PACE.nc', new_nc);

% Add TSS_pred variable with correct _FillValue
nccreate(new_nc, 'TSS_pred', ...
    'Dimensions', {'x', x_dim, 'y', y_dim}, ...
    'Datatype', 'single', ...
    'FillValue', -32767);  % Numeric _FillValue for single

ncwrite(new_nc, 'TSS_pred', tss_map);

% Add attributes to TSS_pred
ncwriteatt(new_nc, 'TSS_pred', 'long_name', 'Predicted Total Suspended Solids');
ncwriteatt(new_nc, 'TSS_pred', 'units', 'mg/L');
ncwriteatt(new_nc, 'TSS_pred', 'coordinates', 'lat lon');
% ncwriteatt(new_nc, 'TSS_pred', '_FillValue', -32767);

% % Update global history
% history = ncreadatt('subset_PACE.nc', '/', 'history');
% new_history = [history; "Predicted TSS added using RF model on " + datestr(now)];
% ncwriteatt(new_nc, '/', 'history', new_history);

% Step 11: Read predicted_tss.nc and visualize TSS
tss_pred_loaded = ncread('Simulated_predicted_tss.nc', 'TSS_pred');
lon_loaded = ncread('Simulated_predicted_tss.nc', 'longitude');
lat_loaded = ncread('Simulated_predicted_tss.nc', 'latitude');

% Bigger bounding box
lon_min = min(lon_loaded(:)) - 0.5;
lon_max = max(lon_loaded(:)) + 0.5;
lat_min = min(lat_loaded(:)) - 0.5;
lat_max = max(lat_loaded(:)) + 0.5;

% Plot using m_map with high-res coastline
m_proj('lambert', 'long', [lon_min lon_max], 'lat', [lat_min lat_max]);
figure(2);
m_pcolor(lon_loaded, lat_loaded, tss_pred_loaded);
% shading flat;
hold on;
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');  % High-resolution coastline and internal waters
hold off;
h = colorbar;
h.Label.String = 'TSS (mg/L)';
h.Label.FontSize = 12;
caxis([0 5]);  % Adjusted range
m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
title('Loaded Predicted TSS (mg/L) from NC File', 'FontSize', 14);
colormap(jet);
saveas(gcf, 'Simulated_tss_map_loaded.png');

disp('Predicted TSS map saved as tss_map.png using m_map Lambert projection with high-res coastline');
disp('New NetCDF file saved as predicted_tss.nc with TSS_pred variable');
disp('Loaded TSS map saved as tss_map_loaded.png using m_map Lambert projection with high-res coastline');
disp('Sample Rrs spectra saved as sample_rrs_spectra.png');

%%

clear; close all; clc;

% Step 1: Load the saved model
load('MODIS_rf_model_best.mat', 'rf_model_best', 'wavelengths');

% Model wavelengths (403 to 718 nm, 117 bands)
model_wl = wavelengths;

% Step 2: Inspect NetCDF file to get Rrs variables and their wavelengths
info = ncinfo('subset_PACE.nc');
var_names = {info.Variables.Name};

% Find Rrs variables (starting with 'Rrs_', excluding 'Rrs_unc_')
rrs_vars = var_names(startsWith(var_names, 'Rrs_') & ~startsWith(var_names, 'Rrs_unc_'));
num_rrs = length(rrs_vars);

% Extract actual wavelengths from attributes
nc_wl = zeros(1, num_rrs);
for i = 1:num_rrs
    nc_wl(i) = ncreadatt('subset_PACE.nc', rrs_vars{i}, 'radiation_wavelength');
end

% Sort by wavelength
[sorted_nc_wl, sort_idx] = sort(nc_wl);
sorted_rrs_vars = rrs_vars(sort_idx);

% Step 3: Map model wavelengths to nearest NC wavelengths/variables
mapped_vars = cell(1, length(model_wl));
mapped_nc_wl = zeros(1, length(model_wl));
for i = 1:length(model_wl)
    [~, idx] = min(abs(sorted_nc_wl - model_wl(i)));
    mapped_vars{i} = sorted_rrs_vars{idx};
    mapped_nc_wl(i) = sorted_nc_wl(idx);
    % fprintf('Model wl %d mapped to NC wl %.3f (%s)\n', model_wl(i), mapped_nc_wl(i), mapped_vars{i});
end

% Identify indices for 650-719 nm bands
high_impact_idx = find(mapped_nc_wl >= 650 & mapped_nc_wl <= 719);
fprintf('High-impact bands (650-719 nm): %d bands\n', length(high_impact_idx));

% Step 4: Read dimensions and lat/lon
x_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'x')).Length;  % 273
y_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'y')).Length;  % 426
lon = ncread('subset_PACE.nc', 'longitude');
lat = ncread('subset_PACE.nc', 'latitude');

% Compute bounding box with padding
lon_min = min(lon(:)) - 0.5;  % Pad by 0.5 degrees
lon_max = max(lon(:)) + 0.5;
lat_min = min(lat(:)) - 0.5;
lat_max = max(lat(:)) + 0.5;

% Step 5: Read raw Rrs data into 3D array (x, y, wl)
rrs_3d = zeros(x_dim, y_dim, length(mapped_vars), 'double');
for i = 1:length(mapped_vars)
    var = mapped_vars{i};
    rrs_raw = ncread('subset_PACE.nc', var);
    fill = ncreadatt('subset_PACE.nc', var, '_FillValue');
    rrs_raw(rrs_raw == fill) = NaN;
    rrs_3d(:,:,i) = double(rrs_raw);  % Use raw values directly
    % fprintf('Band %s (wl %.3f): raw min %.6f, max %.6f, std %.6f\n', ...
    %     var, mapped_nc_wl(i), min(rrs_raw(:),'omitnan'), max(rrs_raw(:),'omitnan'), ...
    %     std(rrs_raw(:),'omitnan'));
end

% Create mask for NaNs
nan_mask = any(isnan(rrs_3d), 3);

% Debug: Print Rrs for sample pixels
good_pixels = find(~nan_mask);
num_good = length(good_pixels);
fprintf('Sample Rrs for first 3 good pixels (wl 402.654, 555.044, 716.817):\n');
for g = 1:min(3, num_good)
    [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
    % fprintf('Pixel (%d,%d): Rrs_403=%.6f, Rrs_555=%.6f, Rrs_717=%.6f\n', ...
    %     ix, iy, rrs_3d(ix,iy,1), rrs_3d(ix,iy,find(mapped_nc_wl >= 555,1)), rrs_3d(ix,iy,end));
end

% Step 6: Optimization - Identify good pixels
if num_good == 0
    error('No valid pixels found; all masked as invalid.');
end
fprintf('Computing for %d valid pixels (%.1f%% of total).\n', num_good, 100*num_good/(x_dim*y_dim));

% Extract Rrs for good pixels
rrs_2d_good = zeros(num_good, length(model_wl), 'double');
for g = 1:num_good
    [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
    rrs_2d_good(g, :) = squeeze(rrs_3d(ix, iy, :))';
end

% Debug: Check Rrs variance and plot spectra
rrs_std_pixels = std(rrs_2d_good, 1, 1);
rrs_std_wl = std(rrs_2d_good, 1, 2);
fprintf('Rrs std across pixels: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_pixels), max(rrs_std_pixels), mean(rrs_std_pixels));
fprintf('Rrs std across wl: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_wl), max(rrs_std_wl), mean(rrs_std_wl));
figure(3);
hold on;
colors = {'b-', 'r-', 'g-'};
for g = 1:min(3, num_good)
    plot(model_wl, rrs_2d_good(g,:), colors{g}, 'LineWidth', 2, 'DisplayName', sprintf('Pixel %d', g));
end
xlabel('Wavelength (nm)');
ylabel('Rrs (raw units)');
title('Sample Rrs Spectra (First 3 Good Pixels)');
legend('show');
grid on;
saveas(gcf, 'sample_rrs_spectra.png');

% Step 7: Test different weights for 650-719 nm bands
weight_factors = [15, 35, 55];
for w = 1:length(weight_factors)
    weight = weight_factors(w);
    fprintf('\nTesting weight = %d for bands 650-719 nm\n', weight);
    
    % Vectorized second derivatives
    second_deriv_good = zeros(num_good, length(model_wl), 'double');
    smoothed_spectra = zeros(size(rrs_2d_good));
    for g = 1:num_good
        spec = rrs_2d_good(g, :);
        spec = max(spec, eps);
        try
            smoothed_spectra(g, :) = sgolayfilt(spec, 3, 11);
        catch
            smoothed_spectra(g, :) = smoothdata(spec, 'movmean', 5);
        end
    end
    % Debug: Check dimensions before gradient
    fprintf('Size of smoothed_spectra: %s\n', mat2str(size(smoothed_spectra)));
    fprintf('Size of model_wl: %s\n', mat2str(size(model_wl)));
    % Compute gradients along wavelength dimension (second dimension)
    first_deriv = gradient(smoothed_spectra, model_wl, 2);
    second_deriv_good = gradient(first_deriv, model_wl, 2);
    second_deriv_good(:, high_impact_idx) = second_deriv_good(:, high_impact_idx) * weight;
    
    % Debug: Check variance in second derivatives
    deriv_std_pixels = std(second_deriv_good, 1, 1);
    deriv_std_wl = std(second_deriv_good, 1, 2);
    fprintf('Second deriv std across pixels (weight=%d): min %.6f, max %.6f, mean %.6f\n', ...
        weight, min(deriv_std_pixels), max(deriv_std_pixels), mean(deriv_std_pixels));
    fprintf('Second deriv std across wl (weight=%d): min %.6f, max %.6f, mean %.6f\n', ...
        weight, min(deriv_std_wl), max(deriv_std_wl), mean(deriv_std_wl));
    fprintf('Sample second deriv for first pixel (weight=%d): min %.6f, max %.6f\n', ...
        weight, min(second_deriv_good(1,:)), max(second_deriv_good(1,:)));
    
    % Predict TSS
    tss_pred_good = predict(rf_model_best, second_deriv_good);
    
    % Use predictions directly (no distance-based scaling)
    tss_pred_scaled = tss_pred_good;
    
    % Debug: Check predictions
    fprintf('Predicted TSS (weight=%d) before clip min: %.4f, max: %.4f, mean: %.4f, std: %.4f\n', ...
        weight, min(tss_pred_scaled), max(tss_pred_scaled), mean(tss_pred_scaled), std(tss_pred_scaled));
    
    % Create TSS map
    tss_map = NaN(x_dim, y_dim);
    tss_map(good_pixels) = tss_pred_scaled;
    tss_map(~isnan(tss_map)) = max(0, min(100, tss_map(~isnan(tss_map))));
    
    % Plot TSS map
    figure(10 + w);
    m_proj('lambert', 'long', [lon_min lon_max], 'lat', [lat_min lat_max]);
    m_pcolor(lon, lat, tss_map);
    shading flat;
    hold on;
    m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');
    hold off;
    h = colorbar;
    h.Label.String = 'TSS (mg/L)';
    h.Label.FontSize = 12;
    caxis([0 50]);
    m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
    title(sprintf('Predicted TSS (mg/L), Weight=%d for 650-719 nm', weight), 'FontSize', 14);
    colormap(jet);
    saveas(gcf, sprintf('tss_map_weight_%d.png', weight));
end

% Step 8: Save final TSS (weight=5) to NetCDF
tss_map = NaN(x_dim, y_dim);
tss_map(good_pixels) = tss_pred_scaled;  % Use last weight (5)
tss_map(~isnan(tss_map)) = max(0, min(100, tss_map(~isnan(tss_map))));
new_nc = 'predicted_tss.nc';
if exist(new_nc, 'file')
    delete(new_nc);
end
nccreate(new_nc, 'longitude', 'Dimensions', {'x', x_dim, 'y', y_dim}, 'Datatype', 'double');
nccreate(new_nc, 'latitude', 'Dimensions', {'x', x_dim, 'y', y_dim}, 'Datatype', 'double');
nccreate(new_nc, 'TSS_pred', 'Dimensions', {'x', x_dim, 'y', y_dim}, 'Datatype', 'single', 'FillValue', -32767);
ncwrite(new_nc, 'longitude', lon);
ncwrite(new_nc, 'latitude', lat);
ncwrite(new_nc, 'TSS_pred', tss_map);
ncwriteatt(new_nc, 'TSS_pred', 'long_name', 'Predicted Total Suspended Solids');
ncwriteatt(new_nc, 'TSS_pred', 'units', 'mg/L');
ncwriteatt(new_nc, 'TSS_pred', 'coordinates', 'latitude longitude');
% ncwriteatt(new_nc, 'TSS_pred', '_FillValue', single(-32767));
ncwriteatt(new_nc, '/', 'history', sprintf('Predicted TSS added using RF model on %s', datestr(now)));

% Step 9: Read and visualize TSS
tss_pred_loaded = ncread(new_nc, 'TSS_pred');
lon_loaded = ncread(new_nc, 'longitude');
lat_loaded = ncread(new_nc, 'latitude');
lon_min = min(lon_loaded(:)) - 0.5;
lon_max = max(lon_loaded(:)) + 0.5;
lat_min = min(lat_loaded(:)) - 0.5;
lat_max = max(lat_loaded(:)) + 0.5;
figure(2);
m_proj('lambert', 'long', [lon_min lon_max], 'lat', [lat_min lat_max]);
m_pcolor(lon_loaded, lat_loaded, tss_pred_loaded);
shading flat;
hold on;
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');
hold off;
h = colorbar;
h.Label.String = 'TSS (mg/L)';
h.Label.FontSize = 12;
caxis([0 50]);
m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
title('Loaded Predicted TSS (mg/L) from NC File', 'FontSize', 14);
colormap(jet);
saveas(gcf, 'tss_map_loaded.png');

disp('Predicted TSS maps saved as tss_map_weight_X.png for weights 1, 2, 5');
disp('New NetCDF file saved as predicted_tss.nc with TSS_pred (weight=5)');
disp('Loaded TSS map saved as tss_map_loaded.png');
disp('Sample Rrs spectra saved as sample_rrs_spectra.png');

%%

clear; close all; clc;

% Step 1: Load the saved model
load('Simulated_PACE_rf_model_raw.mat', 'rf_model_best', 'wavelengths');

% Model wavelengths (403 to 718 nm, 117 bands)
model_wl = wavelengths;

% Step 2: Inspect NetCDF file to get Rrs variables and their wavelengths
info = ncinfo('subset_PACE.nc');
var_names = {info.Variables.Name};

% Find Rrs variables (starting with 'Rrs_', excluding 'Rrs_unc_')
rrs_vars = var_names(startsWith(var_names, 'Rrs_') & ~startsWith(var_names, 'Rrs_unc_'));
num_rrs = length(rrs_vars);

% Extract actual wavelengths from attributes
nc_wl = zeros(1, num_rrs);
for i = 1:num_rrs
    nc_wl(i) = ncreadatt('subset_PACE.nc', rrs_vars{i}, 'radiation_wavelength');
end

% Sort by wavelength
[sorted_nc_wl, sort_idx] = sort(nc_wl);
sorted_rrs_vars = rrs_vars(sort_idx);

% Step 3: Map model wavelengths to nearest NC wavelengths/variables
mapped_vars = cell(1, length(model_wl));
mapped_nc_wl = zeros(1, length(model_wl));
for i = 1:length(model_wl)
    [~, idx] = min(abs(sorted_nc_wl - model_wl(i)));
    mapped_vars{i} = sorted_rrs_vars{idx};
    mapped_nc_wl(i) = sorted_nc_wl(idx);
    % fprintf('Model wl %d mapped to NC wl %.3f (%s)\n', model_wl(i), mapped_nc_wl(i), mapped_vars{i});
end

% Identify indices for 650-719 nm bands
high_impact_idx = find(mapped_nc_wl >= 650 & mapped_nc_wl <= 719);
fprintf('High-impact bands (650-719 nm): %d bands\n', length(high_impact_idx));

% Step 4: Read dimensions and lat/lon
x_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'x')).Length;  % 273
y_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'y')).Length;  % 426
lon = ncread('subset_PACE.nc', 'longitude');
lat = ncread('subset_PACE.nc', 'latitude');

% Compute bounding box with padding
lon_min = min(lon(:)) - 0.5;  % Pad by 0.5 degrees
lon_max = max(lon(:)) + 0.5;
lat_min = min(lat(:)) - 0.5;
lat_max = max(lat(:)) + 0.5;

% Step 5: Read raw Rrs data into 3D array (x, y, wl)
rrs_3d = zeros(x_dim, y_dim, length(mapped_vars), 'double');
for i = 1:length(mapped_vars)
    var = mapped_vars{i};
    rrs_raw = ncread('subset_PACE.nc', var);
    fill = ncreadatt('subset_PACE.nc', var, '_FillValue');
    rrs_raw(rrs_raw == fill) = NaN;
    rrs_3d(:,:,i) = double(rrs_raw);  % Use raw values directly
    % fprintf('Band %s (wl %.3f): raw min %.6f, max %.6f, std %.6f\n', ...
    %     var, mapped_nc_wl(i), min(rrs_raw(:),'omitnan'), max(rrs_raw(:),'omitnan'), ...
    %     std(rrs_raw(:),'omitnan'));
end

% Create mask for NaNs
nan_mask = any(isnan(rrs_3d), 3);

% Debug: Print Rrs for sample pixels
good_pixels = find(~nan_mask);
num_good = length(good_pixels);
fprintf('Sample Rrs for first 3 good pixels (wl 402.654, 555.044, 716.817):\n');
for g = 1:min(3, num_good)
    [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
    % fprintf('Pixel (%d,%d): Rrs_403=%.6f, Rrs_555=%.6f, Rrs_717=%.6f\n', ...
    %     ix, iy, rrs_3d(ix,iy,1), rrs_3d(ix,iy,find(mapped_nc_wl >= 555,1)), rrs_3d(ix,iy,end));
end

% Step 6: Optimization - Identify good pixels
if num_good == 0
    error('No valid pixels found; all masked as invalid.');
end
fprintf('Computing for %d valid pixels (%.1f%% of total).\n', num_good, 100*num_good/(x_dim*y_dim));

% Extract Rrs for good pixels
rrs_2d_good = zeros(num_good, length(model_wl), 'double');
for g = 1:num_good
    [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
    rrs_2d_good(g, :) = squeeze(rrs_3d(ix, iy, :))';
end

% Debug: Check Rrs variance and plot spectra
rrs_std_pixels = std(rrs_2d_good, 1, 1);
rrs_std_wl = std(rrs_2d_good, 1, 2);
fprintf('Rrs std across pixels: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_pixels), max(rrs_std_pixels), mean(rrs_std_pixels));
fprintf('Rrs std across wl: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_wl), max(rrs_std_wl), mean(rrs_std_wl));
figure(3);
hold on;
colors = {'b-', 'r-', 'g-'};
for g = 1:min(3, num_good)
    plot(model_wl, rrs_2d_good(g,:), colors{g}, 'LineWidth', 2, 'DisplayName', sprintf('Pixel %d', g));
end
xlabel('Wavelength (nm)');
ylabel('Rrs (raw units)');
title('Sample Rrs Spectra (First 3 Good Pixels)');
legend('show');
grid on;
saveas(gcf, 'sample_rrs_spectra.png');

% Step 7: Test different weights for 650-719 nm bands
weight_factors = [1, 2, 5];
for w = 1:length(weight_factors)
    weight = weight_factors(w);
    fprintf('\nTesting weight = %d for bands 650-719 nm\n', weight);
    
    % Vectorized second derivatives
    second_deriv_good = zeros(num_good, length(model_wl), 'double');
    smoothed_spectra = zeros(size(rrs_2d_good));
    for g = 1:num_good
        spec = rrs_2d_good(g, :);
        spec = max(spec, eps);
        try
            smoothed_spectra(g, :) = sgolayfilt(spec, 3, 11);
        catch
            smoothed_spectra(g, :) = smoothdata(spec, 'movmean', 5);
        end
    end
    % Debug: Check dimensions before gradient
    fprintf('Size of smoothed_spectra: %s\n', mat2str(size(smoothed_spectra)));
    fprintf('Size of model_wl: %s\n', mat2str(size(model_wl)));
    % Compute gradients along wavelength dimension (second dimension)
    first_deriv = gradient(smoothed_spectra, model_wl, 2);
    second_deriv_good = gradient(first_deriv, model_wl, 2);
    second_deriv_good(:, high_impact_idx) = second_deriv_good(:, high_impact_idx) * weight;
    
    % Debug: Check variance in second derivatives
    deriv_std_pixels = std(second_deriv_good, 1, 1);
    deriv_std_wl = std(second_deriv_good, 1, 2);
    fprintf('Second deriv std across pixels (weight=%d): min %.6f, max %.6f, mean %.6f\n', ...
        weight, min(deriv_std_pixels), max(deriv_std_pixels), mean(deriv_std_pixels));
    fprintf('Second deriv std across wl (weight=%d): min %.6f, max %.6f, mean %.6f\n', ...
        weight, min(deriv_std_wl), max(deriv_std_wl), mean(deriv_std_wl));
    fprintf('Sample second deriv for first pixel (weight=%d): min %.6f, max %.6f\n', ...
        weight, min(second_deriv_good(1,:)), max(second_deriv_good(1,:)));
    
    % Predict TSS
    tss_pred_good = predict(rf_model_best, second_deriv_good);
    
    % Use predictions directly (no distance-based scaling)
    tss_pred_scaled = tss_pred_good;
    
    % Debug: Check predictions
    fprintf('Predicted TSS (weight=%d) before clip min: %.4f, max: %.4f, mean: %.4f, std: %.4f\n', ...
        weight, min(tss_pred_scaled), max(tss_pred_scaled), mean(tss_pred_scaled), std(tss_pred_scaled));
    
    % Create TSS map
    tss_map = NaN(x_dim, y_dim);
    tss_map(good_pixels) = tss_pred_scaled;
    tss_map(~isnan(tss_map)) = max(0, min(100, tss_map(~isnan(tss_map))));
    
    % Plot TSS map
    figure(10 + w);
    m_proj('lambert', 'long', [lon_min lon_max], 'lat', [lat_min lat_max]);
    m_pcolor(lon, lat, tss_map);
    shading flat;
    hold on;
    m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');
    hold off;
    h = colorbar;
    h.Label.String = 'TSS (mg/L)';
    h.Label.FontSize = 12;
    caxis([0 50]);
    m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
    title(sprintf('Predicted TSS (mg/L), Weight=%d for 650-719 nm', weight), 'FontSize', 14);
    colormap(jet);
    saveas(gcf, sprintf('tss_map_weight_%d.png', weight));
end

% Step 8: Save final TSS (weight=5) to NetCDF
tss_map = NaN(x_dim, y_dim);
tss_map(good_pixels) = tss_pred_scaled;  % Use last weight (5)
tss_map(~isnan(tss_map)) = max(0, min(100, tss_map(~isnan(tss_map))));
new_nc = 'predicted_tss.nc';
if exist(new_nc, 'file')
    delete(new_nc);
end
nccreate(new_nc, 'longitude', 'Dimensions', {'x', x_dim, 'y', y_dim}, 'Datatype', 'double');
nccreate(new_nc, 'latitude', 'Dimensions', {'x', x_dim, 'y', y_dim}, 'Datatype', 'double');
nccreate(new_nc, 'TSS_pred', 'Dimensions', {'x', x_dim, 'y', y_dim}, 'Datatype', 'single', 'FillValue', -32767);
ncwrite(new_nc, 'longitude', lon);
ncwrite(new_nc, 'latitude', lat);
ncwrite(new_nc, 'TSS_pred', tss_map);
ncwriteatt(new_nc, 'TSS_pred', 'long_name', 'Predicted Total Suspended Solids');
ncwriteatt(new_nc, 'TSS_pred', 'units', 'mg/L');
ncwriteatt(new_nc, 'TSS_pred', 'coordinates', 'latitude longitude');
%vncwriteatt(new_nc, 'TSS_pred', '_FillValue', single(-32767));
ncwriteatt(new_nc, '/', 'history', sprintf('Predicted TSS added using RF model on %s', datestr(now)));

% Step 9: Read and visualize TSS
tss_pred_loaded = ncread(new_nc, 'TSS_pred');
lon_loaded = ncread(new_nc, 'longitude');
lat_loaded = ncread(new_nc, 'latitude');
lon_min = min(lon_loaded(:)) - 0.5;
lon_max = max(lon_loaded(:)) + 0.5;
lat_min = min(lat_loaded(:)) - 0.5;
lat_max = max(lat_loaded(:)) + 0.5;

figure(2);
m_proj('lambert', 'long', [lon_min lon_max], 'lat', [lat_min lat_max]);
m_pcolor(lon_loaded, lat_loaded, tss_pred_loaded);
shading flat;
hold on;
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');
hold off;
h = colorbar;
h.Label.String = 'TSS (mg/L)';
h.Label.FontSize = 12;
caxis([0 50]);
m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
title('Loaded Predicted TSS (mg/L) from NC File', 'FontSize', 14);
colormap(jet);
saveas(gcf, 'tss_map_loaded.png');

disp('Predicted TSS maps saved as tss_map_weight_X.png for weights 1, 2, 5');
disp('New NetCDF file saved as predicted_tss.nc with TSS_pred (weight=5)');
disp('Loaded TSS map saved as tss_map_loaded.png');
disp('Sample Rrs spectra saved as sample_rrs_spectra.png');
