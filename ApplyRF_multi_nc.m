% Apply my model for multiple nc files
% %% 
% clear; close all; clc;
% 
% % Step 1: Load the saved model
% load('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% % load('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', 'rf_model', 'wavelengths'); 
% %load('PACE_rf_2nd_der_insitu20_Sim_validated.mat', 'rf_model', 'wavelengths');
% % Model wavelengths (403 to 718 nm, 117 bands)
% model_wl = wavelengths;
% PACE_l2 = "/Volumes/SRC_HDD_Mas/Data/PACE_Cruise_Matchups/Backups/PACE_20250430T175451.L2_AOP_SAB.nc";
% % PACE_l2 = "/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/yellow_sea/PACE_OCI.20250426T040501.L2.OC_AOP.V3_0_subset.nc";
% % PACE_l2 = "subset_PACE.nc";
% 
% % Step 2: Inspect NetCDF file to get Rrs variables and their wavelengths
% % info = ncinfo('subset_PACE.nc'); 
% info = ncinfo(PACE_l2);
% var_names = {info.Variables.Name};
% 
% % Find Rrs variables (starting with 'Rrs_', excluding 'Rrs_unc_')
% rrs_vars = var_names(startsWith(var_names, 'Rrs_') & ~startsWith(var_names, 'Rrs_unc_'));
% num_rrs = length(rrs_vars);
% 
% % Extract actual wavelengths from attributes
% nc_wl = zeros(1, num_rrs);
% for i = 1:num_rrs
%     nc_wl(i) = ncreadatt(PACE_l2, rrs_vars{i}, 'radiation_wavelength');
% end
% 
% % Sort by wavelength
% [sorted_nc_wl, sort_idx] = sort(nc_wl);
% sorted_rrs_vars = rrs_vars(sort_idx);
% 
% % Step 3: Map model wavelengths to nearest NC wavelengths/variables
% mapped_vars = cell(1, length(model_wl));
% mapped_nc_wl = zeros(1, length(model_wl));
% for i = 1:length(model_wl)
%     [~, idx] = min(abs(sorted_nc_wl - model_wl(i)));
%     mapped_vars{i} = sorted_rrs_vars{idx};
%     mapped_nc_wl(i) = sorted_nc_wl(idx);
%     % fprintf('Model wl %d mapped to NC wl %.3f (%s)\n', model_wl(i), mapped_nc_wl(i), mapped_vars{i});
% end
% 
% % Step 4: Read dimensions and lat/lon
% x_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'x')).Length;
% y_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'y')).Length;
% lon = ncread(PACE_l2, 'longitude');
% lat = ncread(PACE_l2, 'latitude');
% 
% % Compute bounding box with padding
% lon_min = min(lon(:)) - 0.5;  % Pad by 0.5 degrees
% lon_max = max(lon(:)) + 0.5;
% lat_min = min(lat(:)) - 0.5;
% lat_max = max(lat(:)) + 0.5;
% 
% % Step 5: Read raw Rrs data into 3D array (x, y, wl)
% rrs_3d = zeros(x_dim, y_dim, length(mapped_vars), 'double');
% for i = 1:length(mapped_vars)
%     var = mapped_vars{i};
%     rrs_raw = ncread(PACE_l2, var);
%     % Apply fill value
%     fill = ncreadatt(PACE_l2, var, '_FillValue');
%     rrs_raw(rrs_raw == fill) = NaN;
%     % Apply scaling and offset
%     scale = ncreadatt(PACE_l2, var, 'scale_factor');
%     offset = ncreadatt(PACE_l2, var, 'add_offset');
%     rrs_scaled = double(rrs_raw); % No need for scale factor
%     % rrs_scaled = double(rrs_raw) * scale + offset;
%     rrs_3d(:,:,i) = rrs_scaled;
% end
% 
% % Create mask for NaNs
% nan_mask = any(isnan(rrs_3d), 3);
% % % Check for l2_flags to improve masking % don't apply this for my data
% % try
% %     l2_flags = ncread(PACE_l2, 'l2_flags');
% %     invalid_mask = bitand(l2_flags, 1+2+4) > 0;  % Example: land (1), cloud (2), fail (4)
% %     nan_mask = nan_mask | invalid_mask;
% %     disp('Applied l2_flags masking for land/clouds');
% % catch
% %     disp('l2_flags not found; using only Rrs NaN mask');
% % end
% 
% % Debug: Print Rrs for sample pixels
% good_pixels = find(~nan_mask);
% num_good = length(good_pixels);
% fprintf('Sample Rrs for first 3 good pixels (wl 402.654, 555.044, 716.817):\n');
% for g = 1:min(3, num_good)
%     [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
%     % fprintf('Pixel (%d,%d): Rrs_403=%.6f, Rrs_555=%.6f, Rrs_717=%.6f\n', ...
%     %     ix, iy, rrs_3d(ix,iy,1), rrs_3d(ix,iy,find(mapped_nc_wl >= 555,1)), rrs_3d(ix,iy,end));
% end
% 
% 
% % Step 6: Optimization - Identify good pixels
% if num_good == 0
%     error('No valid pixels found; all masked as invalid.');
% end
% fprintf('Computing for %d valid pixels (%.1f%% of total).\n', num_good, 100*num_good/(x_dim*y_dim));
% 
% % Extract Rrs for good pixels
% rrs_2d_good = zeros(num_good, length(model_wl), 'double');
% for g = 1:num_good
%     [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
%     rrs_2d_good(g, :) = squeeze(rrs_3d(ix, iy, :))';
% end
% 
% % Debug: Check Rrs variance and plot spectra
% rrs_std_pixels = std(rrs_2d_good, 1, 1);
% rrs_std_wl = std(rrs_2d_good, 1, 2);
% fprintf('Rrs std across pixels: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_pixels), max(rrs_std_pixels), mean(rrs_std_pixels));
% fprintf('Rrs std across wl: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_wl), max(rrs_std_wl), mean(rrs_std_wl));
% 
% % Step 7: Compute second derivatives for good pixels
% second_deriv_good = zeros(num_good, length(model_wl), 'double');
% smoothed_spectra = zeros(size(rrs_2d_good));
% for g = 1:num_good
%     spec = rrs_2d_good(g, :);
%     spec = max(spec, eps);  % Ensure positive
%     try
%         smoothed_spectra(g, :) = sgolayfilt(spec, 3, 11);  % 3rd-order, 11-point window
%     catch
%         smoothed_spectra(g, :) = smoothdata(spec, 'movmean', 5);  % Fallback
%     end
% end
% % Compute gradients along wavelength dimension
% first_deriv = gradient(smoothed_spectra, model_wl, 2);
% second_deriv_good = gradient(first_deriv, model_wl, 2);
% 
% % Debug: Check variance in second derivatives
% deriv_std_pixels = std(second_deriv_good, 1, 1);
% deriv_std_wl = std(second_deriv_good, 1, 2);
% fprintf('Second deriv std across pixels: min %.6f, max %.6f, mean %.6f\n', ...
%     min(deriv_std_pixels), max(deriv_std_pixels), mean(deriv_std_pixels));
% fprintf('Second deriv std across wl: min %.6f, max %.6f, mean %.6f\n', ...
%     min(deriv_std_wl), max(deriv_std_wl), mean(deriv_std_wl));
% fprintf('Sample second deriv for first pixel: min %.6f, max %.6f\n', ...
%     min(second_deriv_good(1,:)), max(second_deriv_good(1,:)));
% 
% % Step 8: Predict TSS for good pixels
% tss_pred_good = predict(rf_model, second_deriv_good);
% 
% % Debug: Check predictions
% fprintf('Predicted TSS before clip min: %.4f, max: %.4f, mean: %.4f, std: %.4f\n', ...
%     min(tss_pred_good), max(tss_pred_good), mean(tss_pred_good), std(tss_pred_good));
% 
% % Create TSS map
% tss_map = NaN(x_dim, y_dim);
% tss_map(good_pixels) = tss_pred_good;
% tss_map(~isnan(tss_map)) = max(0, min(1000, tss_map(~isnan(tss_map))));  % Clip to [0, 300]
% 
% %% Plot TSS map
% 
% figure(10);
% m_proj('lambert', 'long', [-81.7 -80.], 'lat', [30.3 33.2]);
% % m_proj('lambert', 'long', [-51 -46], 'lat', [-2 4.5]);% Amazon
% % m_proj('lambert', 'long', [86 92], 'lat', [18 23]);% BoB
% %m_proj('lambert', 'long', [116 127], 'lat', [32 42]);% Yellow Sea
% % m_proj('lambert', 'long', [120 124], 'lat', [32 36]);% Yellow Sea
% m_pcolor(lon, lat, tss_map);
% shading flat;
% hold on;
% % m_gshhs_h('patch', [.7 .7 .7], 'edgecolor', 'k');
% 
% % High-resolution coastlines (using different GSHHS levels)
% m_gshhs('h', 'patch', [0.85 0.85 0.85], 'edgecolor', 'k', 'linewidth', 1.5); % High-res coast
% % Add filename text in a nice box
% filename = 'PACE_20250430T175451.L2_AOP_SAB.nc'; % Your filename here
% m_text(-81.5, 31.6, 'OCI.20250306T181826' , ...
%        'FontSize', 8, 'FontWeight', 'normal', ...
%        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
%        'Margin', 1, 'HorizontalAlignment', 'left');
% hold off;
% h = colorbar;
% h.Label.String = 'TSS (mg L^{-1})';
% h.Label.FontSize = 12;
% m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
% title('Predicted TSS (mg L^{-1}) Model', 'FontSize', 14);
% 
% caxis([0 30]);
% m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
% title('Predicted TSS (mg L^{-1}) ', 'FontSize', 14);
% 
% colormap(jet);
% figs_name = info.Filename(end-39:end);
% saveas(gcf, '/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/Figures',  figs_name + 'Simulated_Altamaha.png');

%% For multiple files
% Apply my model for multiple nc files
%% 
clear; clc; close all% No close all here, we'll handle figures inside the loop

% Step 1: Load the saved model (outside loop, same for all files)
% load('PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_noSmoothing.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% load('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
load('/Users/masud/OneDriveUGA/QWIP/New_350_719/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');

% load('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', 'rf_model', 'wavelengths'); 
% load('PACE_rf_2nd_der_insitu20_Sim_validated.mat', 'rf_model', 'wavelengths');
% Model wavelengths (403 to 718 nm, 117 bands)
model_wl = wavelengths; 

% Define input and output directories
input_dir = '/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional'; % Adjust if needed to include multiple subdirs or use a different path
output_dir = '/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional';
% Delete all temporary files
%find /Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/data/Altamaha_SAB -type f -name "._*" -delete % Terminal code
% Get list of .nc files (you can modify to include multiple directories or specific files)
nc_files = dir(fullfile(input_dir, '*.nc'));
% Alternatively, hardcode a list if specific files are needed:
% nc_files = {...
%     fullfile('/Volumes/SRC_HDD_Mas/Data/PACE_Cruise_Matchups/Backups/', 'PACE_20250430T175451.L2_AOP_SAB.nc'), ...
%     fullfile('/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/yellow_sea/', 'PACE_OCI.20250426T040501.L2.OC_AOP.V3_0_subset.nc'), ...
%     'subset_PACE.nc' ...
% };
% Then in loop: PACE_l2 = nc_files{f};

%% Loop over each .nc file
for f =  1 :length(nc_files)
    PACE_l2 = fullfile(input_dir, nc_files(f).name);
    % PACE_l2 = "subset_PACE.nc";
    [~, basename, ~] = fileparts(PACE_l2); % Extract basename for saving and text
    
    fprintf('Processing file: %s\n', PACE_l2);
    
    % Step 2: Inspect NetCDF file to get Rrs variables and their wavelengths
    info = ncinfo(PACE_l2);
    var_names = {info.Variables.Name};

    % Find Rrs variables (starting with 'Rrs_', excluding 'Rrs_unc_')
    rrs_vars = var_names(startsWith(var_names, 'Rrs_') & ~startsWith(var_names, 'Rrs_unc_'));
    num_rrs = length(rrs_vars);

    % Extract actual wavelengths from attributes
    nc_wl = zeros(1, num_rrs);
    for i = 1:num_rrs
        nc_wl(i) = ncreadatt(PACE_l2, rrs_vars{i}, 'radiation_wavelength');
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

    % Step 4: Read dimensions and lat/lon
    x_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'x')).Length;
    y_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'y')).Length;
    lon = ncread(PACE_l2, 'longitude');
    lat = ncread(PACE_l2, 'latitude');

    % Compute bounding box with padding based on min/max from this file
    lon_min = min(lon(:)) ;  % Pad by 0.5 degrees
    lon_max = max(lon(:)) ;
    lat_min = min(lat(:)) ;
    lat_max = max(lat(:));

    % Step 5: Read raw Rrs data into 3D array (x, y, wl)
    rrs_3d = zeros(x_dim, y_dim, length(mapped_vars), 'double');
    for i = 1:length(mapped_vars)
        var = mapped_vars{i};
        rrs_raw = ncread(PACE_l2, var);
        % Apply fill value
        fill = ncreadatt(PACE_l2, var, '_FillValue');
        rrs_raw(rrs_raw == fill) = NaN;
        % Apply scaling and offset % No need for this case
%         scale = ncreadatt(PACE_l2, var, 'scale_factor');
%         offset = ncreadatt(PACE_l2, var, 'add_offset');
        rrs_scaled = double(rrs_raw); % No need for scale factor as per original
        % rrs_scaled = double(rrs_raw) * scale + offset;
        rrs_3d(:,:,i) = rrs_scaled;
    end

    % Create mask for NaNs
    nan_mask = any(isnan(rrs_3d), 3);
    % % Check for l2_flags to improve masking % don't apply this for my data
%     try
%         l2_flags = ncread(PACE_l2, 'l2_flags');
%         invalid_mask = bitand(l2_flags, 2) > 0;  % Example: land (1), cloud (2), fail (4)
%         nan_mask = nan_mask | invalid_mask;
%         disp('Applied l2_flags masking for land/clouds');
%     catch
%         disp('l2_flags not found; using only Rrs NaN mask');
%     end

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
        warning('No valid pixels found in %s; skipping.', PACE_l2);
        continue;
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

    % Step 7: Compute second derivatives for good pixels
    second_deriv_good = zeros(num_good, length(model_wl), 'double');
    smoothed_spectra = zeros(size(rrs_2d_good));
    for g = 1:num_good
        spec = rrs_2d_good(g, :);
        spec = max(spec, eps);  % Ensure positive
        smoothed_spectra(g, :) = sgolayfilt(spec, 3, 21);
%         smoothed_spectra(g, :) = smoothdata(spec, 'movmean', 5);  % Fallback
%         try
%             smoothed_spectra(g, :) = sgolayfilt(spec, 3, 15);  % 3rd-order, 11-point window
%         catch
%             smoothed_spectra(g, :) = smoothdata(spec, 'movmean', 5);  % Fallback
%         end
    end
    % Compute gradients along wavelength dimension
    first_deriv = gradient(smoothed_spectra, model_wl, 2);
    % first_deriv = gradient(rrs_2d_good, model_wl, 2);
    second_deriv_good = gradient(first_deriv, model_wl, 2);

    % Debug: Check variance in second derivatives
    deriv_std_pixels = std(second_deriv_good, 1, 1);
    deriv_std_wl = std(second_deriv_good, 1, 2);
    fprintf('Second deriv std across pixels: min %.6f, max %.6f, mean %.6f\n', ...
        min(deriv_std_pixels), max(deriv_std_pixels), mean(deriv_std_pixels));
    fprintf('Second deriv std across wl: min %.6f, max %.6f, mean %.6f\n', ...
        min(deriv_std_wl), max(deriv_std_wl), mean(deriv_std_wl));
    fprintf('Sample second deriv for first pixel: min %.6f, max %.6f\n', ...
        min(second_deriv_good(1,:)), max(second_deriv_good(1,:)));

    % Step 8: Predict TSS for good pixels
    tss_pred_good = predict(rf_model, second_deriv_good);

    % Debug: Check predictions
    fprintf('Predicted TSS before clip min: %.4f, max: %.4f, mean: %.4f, std: %.4f\n', ...
        min(tss_pred_good), max(tss_pred_good), mean(tss_pred_good), std(tss_pred_good));

    % Create TSS map
    tss_map = NaN(x_dim, y_dim);
    tss_map(good_pixels) = tss_pred_good;
    tss_map(~isnan(tss_map)) = max(0, min(1000, tss_map(~isnan(tss_map))));  % Clip to [0, 1000]
    % Plot TSS map with logarithmic color scale, smooth colors, no patch outlines, and NaNs preserved as transparent
% figure(); % New figure for each file
% 
% % Set up Lambert projection for South Atlantic Bight (SAB)
% m_proj('lambert', 'long', [-81.7 -80.], 'lat', [30.3 32]);
% 
% % Create a mask for valid (non-NaN) TSS values
% valid_mask = ~isnan(tss_map);
% 
% % Ensure only positive valid data for logarithmic scale
% valid_tss = tss_map(valid_mask);
% if isempty(valid_tss) || all(valid_tss <= 0)
%     warning('No valid positive TSS data for station %s. Skipping plot.', basename);
%     close(gcf);
%     return;
% end
% 
% % Plot the TSS map with logarithmic color scale and no edge outlines
% h_pcolor = m_pcolor(lon, lat, tss_map);
% set(h_pcolor, 'FaceAlpha', 'flat', 'AlphaData', valid_mask, 'EdgeColor', 'none'); % Make NaN pixels transparent, remove edges
% shading interp; % Use interpolated shading for smooth colors
% set(gca, 'ColorScale', 'log'); % Set logarithmic color scale
% hold on;
% 
% % Add high-resolution coastlines
% m_gshhs('h', 'patch', [0.85 0.85 0.85], 'edgecolor', 'k', 'linewidth', 1.5);
% 
% % Add text annotation with basename
% text_str = basename(10:25); % Adjust as needed for datetime extraction
% m_text(-81.5, 31.8, text_str, ...
%        'FontSize', 8, 'FontWeight', 'normal', ...
%        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
%        'Margin', 1, 'HorizontalAlignment', 'left');
% 
% hold off;
% 
% % Add colorbar and update label
% h = colorbar;
% h.Label.String = 'logTSS (mg L^{-1})';
% h.Label.FontSize = 12;
% caxis([0.05 25]);
% set(h,'FontSize', 12,'FontName','Leelawadee')
% %set(h,'ytick',([0.05 0.1 1 2 5]),'yticklabel',["<0.05", "0.1", "1", "2","5"],'tickdir','out')
% set(get(h,'xlabel'),'string','TSS (mg L^{-1})', 'FontSize', 12,'FontName','Leelawadee')
% 
% % Use jet colormap
% colormap(jet);
% 
% % Add grid and title
% m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
% title('Predicted TSS (log_{10} scale, mg L^{-1})', 'FontSize', 14);
% 
% % Save the figure
% save_filename = fullfile(output_dir, [basename '_Simulated_TSS_LogScale.png']);
% saveas(gcf, save_filename);
% fprintf('Saved map to: %s\n', save_filename);

    % 
    %figure();
    figure('Visible', 'off'); % Create figure without displaying it
    m_proj('lambert', 'long', [-82. -78.0], 'lat', [28.0 34.0]); % Full SAB
    % m_proj('lambert', 'long', [-81.7 -80.], 'lat', [30.3 32]); %Altamaha River
    % m_proj('lambert', 'long', [-50.5 -47], 'lat', [-2 2.5]);% Amazon
    % m_proj('lambert', 'long', [86 92], 'lat', [18 23]);% BoB
    % m_proj('lambert', 'long', [94 96], 'lat', [14.5 16.5]);% Irrawaddy river
    %m_proj('lambert', 'long', [116 127], 'lat', [32 42]);% Yellow Sea
    % m_proj('lambert', 'long', [120 124], 'lat', [32 36]);% Yellow Sea
    m_pcolor(lon, lat, tss_map);
    % shading interp;
    shading flat;
    hold on;
    % m_gshhs_h('patch', [.7 .7 .7], 'edgecolor', 'k');
    
    % High-resolution coastlines (using different GSHHS levels)
    m_gshhs('h', 'patch', [0.85 0.85 0.85], 'edgecolor', 'k', 'linewidth', 1.5); % High-res coast
    % Add filename text in a nice box
    % filename = 'PACE_20250430T175451.L2_AOP_SAB.nc'; % Your filename here
    % m_text(-81.5, 31.8, basename(9:25) , ...
    %        'FontSize', 8, 'FontWeight', 'normal', ...
    %        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
    %        'Margin', 1, 'HorizontalAlignment', 'left');
    
    hold off;
    h = colorbar;
    h.Label.String = 'TSS (mg L^{-1})';
    h.Label.FontSize = 12;
    m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
    title('Predicted TSS (mg L^{-1}) Model', 'FontSize', 14);
    text_str = basename(10:25); % Adjust as needed for datetime extraction
    m_text(-81.5, 31.8, text_str, ...
           'FontSize', 8, 'FontWeight', 'normal', ...
           'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
           'Margin', 1, 'HorizontalAlignment', 'left');
    caxis([0 25]);
    m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
    title('Predicted TSS (mg L^{-1}) ', 'FontSize', 14);
    
    colormap(jet);
    figs_name = info.Filename(end-39:end);
    
    % Save the figure
    save_filename = fullfile(output_dir, [basename 'Amazon_Simu_TSS.png']);
    saveas(gcf, save_filename);
    fprintf('Saved map to: %s\n', save_filename);

    % Close figure to free resources
    close(gcf);
    %% Log scale figures for TSS
    % Create figure without displaying it
    figure('Visible', 'off');
    % m_proj('lambert', 'long', [-81.7 -80.], 'lat', [30.3 32]);
    m_proj('lambert', 'long', [-82. -78.0], 'lat', [28.0 34.0]);
    m_pcolor(lon, lat, log10(tss_map)); % Apply log10 transformation to TSS data
    shading flat;
    hold on;
    m_gshhs('h', 'patch', [0.85 0.85 0.85], 'edgecolor', 'k', 'linewidth', 1.5); % High-res coast
    hold off;
    
    % Add colorbar with log scale
    h = colorbar;
    h.Label.String = 'Log_{10} TSS (mg L^{-1})';
    h.Label.FontSize = 12;
    % Set colorbar ticks to correspond to meaningful TSS values (e.g., 0.1, 1, 10, 100 mg/L)
    h.Ticks = log10([0.01 0.1 1 10 100]); % Log-scale ticks
    h.TickLabels = {'0.01','0.1', '1', '10', '100'}; % Readable labels
    caxis(log10([0.01 100])); % Set color axis for log scale (0.1 to 100 mg/L)
    
    % Add grid and title
    m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
    title('Predicted TSS (mg L^{-1}) - Log Scale', 'FontSize', 14);
    text_str = basename(10:25); % Adjust as needed for datetime extraction
    m_text(-81.5, 31.8, text_str, ...
           'FontSize', 8, 'FontWeight', 'normal', ...
           'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
           'Margin', 1, 'HorizontalAlignment', 'left');
    
    % Apply colormap
    colormap(jet);
    
    % Save the figure
    figs_name = info.Filename(end-39:end);
    save_filename = fullfile(output_dir, [basename 'AltamahaSAB_Simu_TSS_logscale.png']);
    saveas(gcf, save_filename);
    fprintf('Saved map to: %s\n', save_filename);
    
    % Close figure to free resources
    close(gcf);
end
%%
% Script to create a movie from all *_Simulated_TSS.png files in a directory
clear; clc;

% Define the directory containing the PNG files
figures_dir = '/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/AmazonRiver/GoodOnes';

% Get list of all PNG files matching the pattern *_Simulated_TSS.png
png_files = dir(fullfile(figures_dir, '*_Simu_TSS.png'));

if isempty(png_files)
    error('No PNG files found in %s', figures_dir);
end

% Extract date-time from filenames for sorting (assuming format like PACE_OCI.YYYYMMDDTHHMMSS...)
% Example: PACE_OCI.20240429T185958.L2.OC_AOP.V3_0_subset_Simulated_TSS.png
% Date string starts after 'PACE_OCI.' and before '.L2'
date_strs = cell(length(png_files), 1);
for i = 1:length(png_files)
    fname = png_files(i).name;
    start_idx = strfind(fname, 'PACE_OCI.') + 9; % After 'PACE_OCI.'
    end_idx = strfind(fname, '.L2') - 1;
    if ~isempty(start_idx) && ~isempty(end_idx)
        date_strs{i} = fname(start_idx:end_idx);
    else
        date_strs{i} = ''; % Fallback if pattern doesn't match
    end
end

% Sort files by date string (chronological order)
[~, sort_idx] = sort(date_strs);
png_files = png_files(sort_idx);

% Define output movie file
output_movie = fullfile(figures_dir, 'TSS_AmazonRiver.mp4');

% Create VideoWriter object (MP4 format, adjust 'FrameRate' as needed)
writerObj = VideoWriter(output_movie, 'MPEG-4');
writerObj.FrameRate = 1; % Frames per second (e.g., 5 for slower animation)
open(writerObj);

% Loop through each PNG file and add to movie
for i = 1:length(png_files)
    fname = fullfile(figures_dir, png_files(i).name);
    fprintf('Adding frame %d/%d: %s\n', i, length(png_files), png_files(i).name);
    
    % Read image
    img = imread(fname);
    
    % Write to video
    writeVideo(writerObj, img);
end

% Close the video writer
close(writerObj);

fprintf('Movie saved to: %s\n', output_movie);
