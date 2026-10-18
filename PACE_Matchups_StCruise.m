%% PACE Validation%
%%  Need to edit this code for PACE Matchups. Copied from : /Users/masud/Documents/GitHub/PACE/ApplyRF_multi_nc.m
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
% % Apply my model for multiple nc files
% %% 
% clear; clc; close all % No close all here, we'll handle figures inside the loop
% 
% % Step 1: Load the saved model (outside loop, same for all files)
% % load('PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_noSmoothing.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% % load('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% % load('/Users/masud/OneDriveUGA/QWIP/New_350_719/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% load('/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat', 'rf_model', 'wavelengths', 'best_rs', 'best_leaf', 'best_trees');
% 
% % load('PACE_rf_model_2nd_deriv_best_Simulated_insitu_validated.mat', 'rf_model', 'wavelengths'); 
% % load('PACE_rf_2nd_der_insitu20_Sim_validated.mat', 'rf_model', 'wavelengths');
% % Model wavelengths (403 to 718 nm, 117 bands)
% model_wl = wavelengths; 
% 
% % Define input and output directories
% input_dir = '/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional'; % Adjust if needed to include multiple subdirs or use a different path
% output_dir = '/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional';
% % Delete all temporary files
% %find /Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/data/Altamaha_SAB -type f -name "._*" -delete % Terminal code
% % Get list of .nc files (you can modify to include multiple directories or specific files)
% nc_files = dir(fullfile(input_dir, '*.nc'));
% % Alternatively, hardcode a list if specific files are needed:
% % nc_files = {...
% %     fullfile('/Volumes/SRC_HDD_Mas/Data/PACE_Cruise_Matchups/Backups/', 'PACE_20250430T175451.L2_AOP_SAB.nc'), ...
% %     fullfile('/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/yellow_sea/', 'PACE_OCI.20250426T040501.L2.OC_AOP.V3_0_subset.nc'), ...
% %     'subset_PACE.nc' ...
% % };
% % Then in loop: PACE_l2 = nc_files{f};
% 
% %% Loop over each .nc file
% for f =  1 :length(nc_files)
%     PACE_l2 = fullfile(input_dir, nc_files(f).name);
%     % PACE_l2 = "subset_PACE.nc";
%     [~, basename, ~] = fileparts(PACE_l2); % Extract basename for saving and text
% 
%     fprintf('Processing file: %s\n', PACE_l2);
% 
%     % Step 2: Inspect NetCDF file to get Rrs variables and their wavelengths
%     info = ncinfo(PACE_l2);
%     var_names = {info.Variables.Name};
% 
%     % Find Rrs variables (starting with 'Rrs_', excluding 'Rrs_unc_')
%     rrs_vars = var_names(startsWith(var_names, 'Rrs_') & ~startsWith(var_names, 'Rrs_unc_'));
%     num_rrs = length(rrs_vars);
% 
%     % Extract actual wavelengths from attributes
%     nc_wl = zeros(1, num_rrs);
%     for i = 1:num_rrs
%         nc_wl(i) = ncreadatt(PACE_l2, rrs_vars{i}, 'radiation_wavelength');
%     end
% 
%     % Sort by wavelength
%     [sorted_nc_wl, sort_idx] = sort(nc_wl);
%     sorted_rrs_vars = rrs_vars(sort_idx);
% 
%     % Step 3: Map model wavelengths to nearest NC wavelengths/variables
%     mapped_vars = cell(1, length(model_wl));
%     mapped_nc_wl = zeros(1, length(model_wl));
%     for i = 1:length(model_wl)
%         [~, idx] = min(abs(sorted_nc_wl - model_wl(i)));
%         mapped_vars{i} = sorted_rrs_vars{idx};
%         mapped_nc_wl(i) = sorted_nc_wl(idx);
%         % fprintf('Model wl %d mapped to NC wl %.3f (%s)\n', model_wl(i), mapped_nc_wl(i), mapped_vars{i});
%     end
% 
%     % Step 4: Read dimensions and lat/lon
%     x_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'x')).Length;
%     y_dim = info.Dimensions(strcmp({info.Dimensions.Name}, 'y')).Length;
%     lon = ncread(PACE_l2, 'longitude');
%     lat = ncread(PACE_l2, 'latitude');
% 
%     % Compute bounding box with padding based on min/max from this file
%     lon_min = min(lon(:)) ;  % Pad by 0.5 degrees
%     lon_max = max(lon(:)) ;
%     lat_min = min(lat(:)) ;
%     lat_max = max(lat(:));
% 
%     % Step 5: Read raw Rrs data into 3D array (x, y, wl)
%     rrs_3d = zeros(x_dim, y_dim, length(mapped_vars), 'double');
%     for i = 1:length(mapped_vars)
%         var = mapped_vars{i};
%         rrs_raw = ncread(PACE_l2, var);
%         % Apply fill value
%         fill = ncreadatt(PACE_l2, var, '_FillValue');
%         rrs_raw(rrs_raw == fill) = NaN;
%         % Apply scaling and offset % No need for this case
% %         scale = ncreadatt(PACE_l2, var, 'scale_factor');
% %         offset = ncreadatt(PACE_l2, var, 'add_offset');
%         rrs_scaled = double(rrs_raw); % No need for scale factor as per original
%         % rrs_scaled = double(rrs_raw) * scale + offset;
%         rrs_3d(:,:,i) = rrs_scaled;
%     end
% 
%     % Create mask for NaNs
%     nan_mask = any(isnan(rrs_3d), 3);
%     % % Check for l2_flags to improve masking % don't apply this for my data
% %     try
% %         l2_flags = ncread(PACE_l2, 'l2_flags');
% %         invalid_mask = bitand(l2_flags, 2) > 0;  % Example: land (1), cloud (2), fail (4)
% %         nan_mask = nan_mask | invalid_mask;
% %         disp('Applied l2_flags masking for land/clouds');
% %     catch
% %         disp('l2_flags not found; using only Rrs NaN mask');
% %     end
% 
%     % Debug: Print Rrs for sample pixels
%     good_pixels = find(~nan_mask);
%     num_good = length(good_pixels);
%     fprintf('Sample Rrs for first 3 good pixels (wl 402.654, 555.044, 716.817):\n');
%     for g = 1:min(3, num_good)
%         [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
%         % fprintf('Pixel (%d,%d): Rrs_403=%.6f, Rrs_555=%.6f, Rrs_717=%.6f\n', ...
%         %     ix, iy, rrs_3d(ix,iy,1), rrs_3d(ix,iy,find(mapped_nc_wl >= 555,1)), rrs_3d(ix,iy,end));
%     end
% 
%     % Step 6: Optimization - Identify good pixels
%     if num_good == 0
%         warning('No valid pixels found in %s; skipping.', PACE_l2);
%         continue;
%     end
%     fprintf('Computing for %d valid pixels (%.1f%% of total).\n', num_good, 100*num_good/(x_dim*y_dim));
% 
%     % Extract Rrs for good pixels
%     rrs_2d_good = zeros(num_good, length(model_wl), 'double');
%     for g = 1:num_good
%         [ix, iy] = ind2sub([x_dim y_dim], good_pixels(g));
%         rrs_2d_good(g, :) = squeeze(rrs_3d(ix, iy, :))';
%     end
% 
%     % Debug: Check Rrs variance and plot spectra
%     rrs_std_pixels = std(rrs_2d_good, 1, 1);
%     rrs_std_wl = std(rrs_2d_good, 1, 2);
%     fprintf('Rrs std across pixels: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_pixels), max(rrs_std_pixels), mean(rrs_std_pixels));
%     fprintf('Rrs std across wl: min %.6f, max %.6f, mean %.6f\n', min(rrs_std_wl), max(rrs_std_wl), mean(rrs_std_wl));
% 
%     % Step 7: Compute second derivatives for good pixels
%     second_deriv_good = zeros(num_good, length(model_wl), 'double');
%     smoothed_spectra = zeros(size(rrs_2d_good));
%     for g = 1:num_good
%         spec = rrs_2d_good(g, :);
%         spec = max(spec, eps);  % Ensure positive
%         smoothed_spectra(g, :) = sgolayfilt(spec, 3, 21);
% %         smoothed_spectra(g, :) = smoothdata(spec, 'movmean', 5);  % Fallback
% %         try
% %             smoothed_spectra(g, :) = sgolayfilt(spec, 3, 15);  % 3rd-order, 11-point window
% %         catch
% %             smoothed_spectra(g, :) = smoothdata(spec, 'movmean', 5);  % Fallback
% %         end
%     end
%     % Compute gradients along wavelength dimension
%     first_deriv = gradient(smoothed_spectra, model_wl, 2);
%     % first_deriv = gradient(rrs_2d_good, model_wl, 2);
%     second_deriv_good = gradient(first_deriv, model_wl, 2);
% 
%     % Debug: Check variance in second derivatives
%     deriv_std_pixels = std(second_deriv_good, 1, 1);
%     deriv_std_wl = std(second_deriv_good, 1, 2);
%     fprintf('Second deriv std across pixels: min %.6f, max %.6f, mean %.6f\n', ...
%         min(deriv_std_pixels), max(deriv_std_pixels), mean(deriv_std_pixels));
%     fprintf('Second deriv std across wl: min %.6f, max %.6f, mean %.6f\n', ...
%         min(deriv_std_wl), max(deriv_std_wl), mean(deriv_std_wl));
%     fprintf('Sample second deriv for first pixel: min %.6f, max %.6f\n', ...
%         min(second_deriv_good(1,:)), max(second_deriv_good(1,:)));
% 
%     % Step 8: Predict TSS for good pixels
%     tss_pred_good = predict(rf_model, second_deriv_good);
% 
%     % Debug: Check predictions
%     fprintf('Predicted TSS before clip min: %.4f, max: %.4f, mean: %.4f, std: %.4f\n', ...
%         min(tss_pred_good), max(tss_pred_good), mean(tss_pred_good), std(tss_pred_good));
% 
%     % Create TSS map
%     tss_map = NaN(x_dim, y_dim);
%     tss_map(good_pixels) = tss_pred_good;
%     tss_map(~isnan(tss_map)) = max(0, min(1000, tss_map(~isnan(tss_map))));  % Clip to [0, 1000]
%     % Plot TSS map with logarithmic color scale, smooth colors, no patch outlines, and NaNs preserved as transparent
% % figure(); % New figure for each file
% % 
% % % Set up Lambert projection for South Atlantic Bight (SAB)
% % m_proj('lambert', 'long', [-81.7 -80.], 'lat', [30.3 32]);
% % 
% % % Create a mask for valid (non-NaN) TSS values
% % valid_mask = ~isnan(tss_map);
% % 
% % % Ensure only positive valid data for logarithmic scale
% % valid_tss = tss_map(valid_mask);
% % if isempty(valid_tss) || all(valid_tss <= 0)
% %     warning('No valid positive TSS data for station %s. Skipping plot.', basename);
% %     close(gcf);
% %     return;
% % end
% % 
% % % Plot the TSS map with logarithmic color scale and no edge outlines
% % h_pcolor = m_pcolor(lon, lat, tss_map);
% % set(h_pcolor, 'FaceAlpha', 'flat', 'AlphaData', valid_mask, 'EdgeColor', 'none'); % Make NaN pixels transparent, remove edges
% % shading interp; % Use interpolated shading for smooth colors
% % set(gca, 'ColorScale', 'log'); % Set logarithmic color scale
% % hold on;
% % 
% % % Add high-resolution coastlines
% % m_gshhs('h', 'patch', [0.85 0.85 0.85], 'edgecolor', 'k', 'linewidth', 1.5);
% % 
% % % Add text annotation with basename
% % text_str = basename(10:25); % Adjust as needed for datetime extraction
% % m_text(-81.5, 31.8, text_str, ...
% %        'FontSize', 8, 'FontWeight', 'normal', ...
% %        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
% %        'Margin', 1, 'HorizontalAlignment', 'left');
% % 
% % hold off;
% % 
% % % Add colorbar and update label
% % h = colorbar;
% % h.Label.String = 'logTSS (mg L^{-1})';
% % h.Label.FontSize = 12;
% % caxis([0.05 25]);
% % set(h,'FontSize', 12,'FontName','Leelawadee')
% % %set(h,'ytick',([0.05 0.1 1 2 5]),'yticklabel',["<0.05", "0.1", "1", "2","5"],'tickdir','out')
% % set(get(h,'xlabel'),'string','TSS (mg L^{-1})', 'FontSize', 12,'FontName','Leelawadee')
% % 
% % % Use jet colormap
% % colormap(jet);
% % 
% % % Add grid and title
% % m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
% % title('Predicted TSS (log_{10} scale, mg L^{-1})', 'FontSize', 14);
% % 
% % % Save the figure
% % save_filename = fullfile(output_dir, [basename '_Simulated_TSS_LogScale.png']);
% % saveas(gcf, save_filename);
% % fprintf('Saved map to: %s\n', save_filename);
% 
%     % 
%     %figure();
%     figure('Visible', 'off'); % Create figure without displaying it
%     m_proj('lambert', 'long', [-82. -78.0], 'lat', [28.0 34.0]); % Full SAB
%     % m_proj('lambert', 'long', [-81.7 -80.], 'lat', [30.3 32]); %Altamaha River
%     % m_proj('lambert', 'long', [-50.5 -47], 'lat', [-2 2.5]);% Amazon
%     % m_proj('lambert', 'long', [86 92], 'lat', [18 23]);% BoB
%     % m_proj('lambert', 'long', [94 96], 'lat', [14.5 16.5]);% Irrawaddy river
%     %m_proj('lambert', 'long', [116 127], 'lat', [32 42]);% Yellow Sea
%     % m_proj('lambert', 'long', [120 124], 'lat', [32 36]);% Yellow Sea
%     m_pcolor(lon, lat, tss_map);
%     % shading interp;
%     shading flat;
%     hold on;
%     % m_gshhs_h('patch', [.7 .7 .7], 'edgecolor', 'k');
% 
%     % High-resolution coastlines (using different GSHHS levels)
%     m_gshhs('h', 'patch', [0.85 0.85 0.85], 'edgecolor', 'k', 'linewidth', 1.5); % High-res coast
%     % Add filename text in a nice box
%     % filename = 'PACE_20250430T175451.L2_AOP_SAB.nc'; % Your filename here
%     % m_text(-81.5, 31.8, basename(9:25) , ...
%     %        'FontSize', 8, 'FontWeight', 'normal', ...
%     %        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
%     %        'Margin', 1, 'HorizontalAlignment', 'left');
% 
%     hold off;
%     h = colorbar;
%     h.Label.String = 'TSS (mg L^{-1})';
%     h.Label.FontSize = 12;
%     m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
%     title('Predicted TSS (mg L^{-1}) Model', 'FontSize', 14);
%     text_str = basename(10:25); % Adjust as needed for datetime extraction
%     m_text(-81.5, 31.8, text_str, ...
%            'FontSize', 8, 'FontWeight', 'normal', ...
%            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
%            'Margin', 1, 'HorizontalAlignment', 'left');
%     caxis([0 25]);
%     m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
%     title('Predicted TSS (mg L^{-1}) ', 'FontSize', 14);
% 
%     colormap(jet);
%     figs_name = info.Filename(end-39:end);
% 
%     % Save the figure
%     save_filename = fullfile(output_dir, [basename 'Amazon_Simu_TSS.png']);
%     saveas(gcf, save_filename);
%     fprintf('Saved map to: %s\n', save_filename);
% 
%     % Close figure to free resources
%     close(gcf);
%     %% Log scale figures for TSS
%     % Create figure without displaying it
%     figure('Visible', 'off');
%     % m_proj('lambert', 'long', [-81.7 -80.], 'lat', [30.3 32]);
%     m_proj('lambert', 'long', [-82. -78.0], 'lat', [28.0 34.0]);
%     m_pcolor(lon, lat, log10(tss_map)); % Apply log10 transformation to TSS data
%     shading flat;
%     hold on;
%     m_gshhs('h', 'patch', [0.85 0.85 0.85], 'edgecolor', 'k', 'linewidth', 1.5); % High-res coast
%     hold off;
% 
%     % Add colorbar with log scale
%     h = colorbar;
%     h.Label.String = 'Log_{10} TSS (mg L^{-1})';
%     h.Label.FontSize = 12;
%     % Set colorbar ticks to correspond to meaningful TSS values (e.g., 0.1, 1, 10, 100 mg/L)
%     h.Ticks = log10([0.01 0.1 1 10 100]); % Log-scale ticks
%     h.TickLabels = {'0.01','0.1', '1', '10', '100'}; % Readable labels
%     caxis(log10([0.01 100])); % Set color axis for log scale (0.1 to 100 mg/L)
% 
%     % Add grid and title
%     m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 12);
%     title('Predicted TSS (mg L^{-1}) - Log Scale', 'FontSize', 14);
%     text_str = basename(10:25); % Adjust as needed for datetime extraction
%     m_text(-81.5, 31.8, text_str, ...
%            'FontSize', 8, 'FontWeight', 'normal', ...
%            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'black', ...
%            'Margin', 1, 'HorizontalAlignment', 'left');
% 
%     % Apply colormap
%     colormap(jet);
% 
%     % Save the figure
%     figs_name = info.Filename(end-39:end);
%     save_filename = fullfile(output_dir, [basename 'AltamahaSAB_Simu_TSS_logscale.png']);
%     saveas(gcf, save_filename);
%     fprintf('Saved map to: %s\n', save_filename);
% 
%     % Close figure to free resources
%     close(gcf);
% end
% %%
% % Script to create a movie from all *_Simulated_TSS.png files in a directory
% clear; clc;
% 
% % Define the directory containing the PNG files
% figures_dir = '/Volumes/SRC_HDD_Mas/Data/PACE_l2_regional/AmazonRiver/GoodOnes';
% 
% % Get list of all PNG files matching the pattern *_Simulated_TSS.png
% png_files = dir(fullfile(figures_dir, '*_Simu_TSS.png'));
% 
% if isempty(png_files)
%     error('No PNG files found in %s', figures_dir);
% end
% 
% % Extract date-time from filenames for sorting (assuming format like PACE_OCI.YYYYMMDDTHHMMSS...)
% % Example: PACE_OCI.20240429T185958.L2.OC_AOP.V3_0_subset_Simulated_TSS.png
% % Date string starts after 'PACE_OCI.' and before '.L2'
% date_strs = cell(length(png_files), 1);
% for i = 1:length(png_files)
%     fname = png_files(i).name;
%     start_idx = strfind(fname, 'PACE_OCI.') + 9; % After 'PACE_OCI.'
%     end_idx = strfind(fname, '.L2') - 1;
%     if ~isempty(start_idx) && ~isempty(end_idx)
%         date_strs{i} = fname(start_idx:end_idx);
%     else
%         date_strs{i} = ''; % Fallback if pattern doesn't match
%     end
% end
% 
% % Sort files by date string (chronological order)
% [~, sort_idx] = sort(date_strs);
% png_files = png_files(sort_idx);
% 
% % Define output movie file
% output_movie = fullfile(figures_dir, 'TSS_AmazonRiver.mp4');
% 
% % Create VideoWriter object (MP4 format, adjust 'FrameRate' as needed)
% writerObj = VideoWriter(output_movie, 'MPEG-4');
% writerObj.FrameRate = 1; % Frames per second (e.g., 5 for slower animation)
% open(writerObj);
% 
% % Loop through each PNG file and add to movie
% for i = 1:length(png_files)
%     fname = fullfile(figures_dir, png_files(i).name);
%     fprintf('Adding frame %d/%d: %s\n', i, length(png_files), png_files(i).name);
% 
%     % Read image
%     img = imread(fname);
% 
%     % Write to video
%     writeVideo(writerObj, img);
% end
% 
% % Close the video writer
% close(writerObj);
% 
% fprintf('Movie saved to: %s\n', output_movie);

%% 15 October 2026
%% Step 1: Load the Saved RF Model
clc; clear; close all;

% Define file paths
mat_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat';
csv_file = 'St_SRC_Optics_cruise_PACE_Apr2025_Oct2026.csv';

% Load the TreeBagger model and wavelengths array
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load and Prepare Matchup Data from CSV
fprintf('Loading matchup data from: %s\n', csv_file);
matchups = readtable(csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Construct the exact column names needed based on the model's wavelengths
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Verify all required columns exist in the CSV before extracting
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    warning('Some model wavelengths are missing from the CSV! Extracting nearest available instead.');
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
        fprintf('Mapped model wl %d to CSV column %s\n', round(model_wl(m)), expected_cols{m});
    end
end

% Extract the Rrs data into an N-by-117 matrix
rrs_matrix = table2array(matchups(:, expected_cols));
fprintf('Successfully extracted %d Rrs columns matching the model.\n', size(rrs_matrix, 2));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

fprintf('Applying Savitzky-Golay smoothing and calculating 2nd derivatives...\n');
for i = 1:N_matchups
    spec = rrs_matrix(i, :);
    spec = max(spec, eps);  % Ensure positive values
    
    % 3rd-order, 21-point window
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21);
end

% Compute gradients along the wavelength dimension
first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Predict TSS using the RF Model
fprintf('Predicting TSS using Random Forest Model...\n');
tss_predicted = predict(rf_model, second_deriv_matrix);

% Clip predictions to bounds [0, 1000] 
tss_predicted = max(0, min(1000, tss_predicted));

%% Step 5: Statistical Analysis
% Clean data by keeping only valid (non-NaN) paired matchups
valid_idx = ~isnan(tss_insitu) & ~isnan(tss_predicted);
obs = tss_insitu(valid_idx);
pred = tss_predicted(valid_idx);

% Calculate N, RMSE, and R-squared
N = length(obs);
if N == 0
    error('No valid matching pairs of in-situ and predicted TSS found.');
end

RMSE = sqrt(mean((pred - obs).^2));
R_mat = corrcoef(obs, pred);
Rsq = R_mat(1,2)^2;

fprintf('\nValidation Stats:\n N = %d\n R^2 = %.4f\n RMSE = %.4f\n', N, Rsq, RMSE);

%% Step 6: Save Predicted vs Insitu Data to CSV
% Create a table with the results and save it
output_csv_name = 'Predicted_vs_Insitu_TSS_Results_v2.csv';
results_table = table(obs, pred, 'VariableNames', {'InSitu_TSS', 'Predicted_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved predicted vs in-situ TSS to: %s\n', output_csv_name);

%% Step 7: Plot Scatter (Predicted vs. In Situ TSS)
fig = figure;
scatter(obs, pred, 60, 'b', 'filled', 'MarkerEdgeColor', 'k');
hold on;

% Add point numbers next to each scatter point
for i = 1:N
    % Adds text slightly offset from the point so it is readable
    text(obs(i), pred(i), sprintf(' %d', i), 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);
end

% Calculate limits to draw a perfect 1:1 line
min_val = min([obs; pred; 0]);
max_val = max([obs; pred]) * 1.1; % Add 10% padding for visual space
plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1.5);

% Plot formatting
grid on;
box on;
axis([min_val max_val min_val max_val]);
xlabel('In Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
title('RF Model Validation: PACE Matchups', 'FontSize', 14);

% Add statistics text box to the top-left of the plot
stats_str = sprintf('N = %d\nR^2 = %.2f\nRMSE = %.2f mg L^{-1}', N, Rsq, RMSE);
text(min_val + (max_val-min_val)*0.05, max_val - (max_val-min_val)*0.1, ...
     stats_str, 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');

%% Step 8: Save the Figure
output_fig_name = 'RF_Model_Validation_Scatter_v2.png';
% Use exportgraphics for a high-quality save (MATLAB R2020a and later)
try
    exportgraphics(fig, output_fig_name, 'Resolution', 300);
catch
    % Fallback for older MATLAB versions
    saveas(fig, output_fig_name);
end
fprintf('Saved scatter plot figure to: %s\n', output_fig_name);

%% log scale
%% Step 1: Load the Saved RF Model
clc; clear; close all;

% Define file paths
mat_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat';
csv_file = 'St_cruise_PACE_26_09_2026.csv';

% Load the TreeBagger model and wavelengths array
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load and Prepare Matchup Data from CSV
fprintf('Loading matchup data from: %s\n', csv_file);
matchups = readtable(csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Construct the exact column names needed based on the model's wavelengths
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Verify and extract columns
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
    end
end
rrs_matrix = table2array(matchups(:, expected_cols));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21);
end
first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Predict TSS using the RF Model
fprintf('Predicting TSS using Random Forest Model...\n');
tss_rf_raw = predict(rf_model, second_deriv_matrix);
tss_rf_raw = max(0, min(1000, tss_rf_raw));

%% Step 5: Post-Hoc Regional Calibration (Using first 7 points)
% Clean data to ensure no NaNs in the first 7 points
valid_idx = ~isnan(tss_insitu) & ~isnan(tss_rf_raw);
obs = tss_insitu(valid_idx);
pred_raw = tss_rf_raw(valid_idx);
N_total = length(obs);

if N_total < 7
    error('Not enough valid data points to calibrate. Need at least 7.');
end

% Isolate the first 7 points for calibration
calib_obs = obs(1:7);
calib_pred = pred_raw(1:7);

% Fit a linear regression model: obs = m * pred + b
% Using polyfit (degree 1)
p = polyfit(calib_pred, calib_obs, 1);
slope = p(1);
intercept = p(2);

fprintf('\nApplied Regional Calibration based on first 7 points:\n');
fprintf('Corrected TSS = (%.4f * RF_TSS) + %.4f\n', slope, intercept);

% Apply the correction to ALL predictions
pred_calibrated = (pred_raw * slope) + intercept;
pred_calibrated = max(0, pred_calibrated); % Prevent negative TSS values

%% Step 6: Statistical Analysis (on Calibrated Data)
% Calculate N, RMSE, and R-squared for the corrected dataset
RMSE = sqrt(mean((pred_calibrated - obs).^2));
R_mat = corrcoef(obs, pred_calibrated);
Rsq = R_mat(1,2)^2;

fprintf('Validation Stats (Calibrated):\n N = %d\n R^2 = %.4f\n RMSE = %.4f\n', N_total, Rsq, RMSE);

%% Step 7: Save Predicted vs Insitu Data to CSV
output_csv_name = 'Predicted_vs_Insitu_TSS_Calibrated.csv';
results_table = table(obs, pred_raw, pred_calibrated, ...
    'VariableNames', {'InSitu_TSS', 'Raw_RF_TSS', 'Calibrated_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved results to: %s\n', output_csv_name);

%% Step 8: Plot Scatter (Calibrated vs. In Situ TSS)
fig = figure;
scatter(obs, pred_calibrated, 60, 'b', 'filled', 'MarkerEdgeColor', 'k');
hold on;

% Highlight the 7 calibration points in a different color (e.g., Red)
scatter(obs(1:7), pred_calibrated(1:7), 60, 'r', 'filled', 'MarkerEdgeColor', 'k');

% Calculate limits
min_val = min([obs; pred_calibrated]);
if min_val <= 0; min_val = 0.01; else; min_val = min_val * 0.8; end
max_val = max([obs; pred_calibrated]) * 1.3; 

% Draw 1:1 line
plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1.5);

% Add point numbers
for i = 1:N_total
    text(obs(i)*1.05, pred_calibrated(i), sprintf('%d', i), 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);
end

set(gca, 'XScale', 'log', 'YScale', 'log'); 
grid on; box on;
axis([min_val max_val min_val max_val]);
xlabel('In Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Calibrated Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
title('Calibrated RF Model Validation', 'FontSize', 14);

% Add statistics
stats_str = sprintf('N = %d\nR^2 = %.2f\nRMSE = %.2f mg L^{-1}', N_total, Rsq, RMSE);
text_x = min_val * 1.2; 
text_y = max_val * 0.6; 
text(text_x, text_y, stats_str, 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');

% Add legend
legend('Validation Points', 'Calibration Points (First 7)', '1:1 Line', 'Location', 'southeast');

%% Step 9: Save the Figure
output_fig_name = 'RF_Model_Validation_Scatter_Calibrated.png';
try
    exportgraphics(fig, output_fig_name, 'Resolution', 300);
catch
    saveas(fig, output_fig_name);
end

%%
%% Step 1: Load the Saved RF Model
clc; clear; close all;

% Define file paths
mat_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat';
csv_file = 'St_cruise_PACE_26_09_2026.csv';

% Load the TreeBagger model and wavelengths array
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load and Prepare Matchup Data from CSV
fprintf('Loading matchup data from: %s\n', csv_file);
matchups = readtable(csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Construct the exact column names needed based on the model's wavelengths
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Verify and extract matching columns
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
    end
end
rrs_matrix = table2array(matchups(:, expected_cols));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21);
end
first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Predict TSS using the RF Model
fprintf('Predicting raw TSS using Random Forest Model...\n');
tss_rf_raw = predict(rf_model, second_deriv_matrix);
tss_rf_raw = max(0, min(1000, tss_rf_raw));

%% Step 5: Post-Hoc Regional Calibration (Using first 7 points)
valid_idx = ~isnan(tss_insitu) & ~isnan(tss_rf_raw);
obs = tss_insitu(valid_idx);
pred_raw = tss_rf_raw(valid_idx);
N_total = length(obs);

if N_total < 7
    error('Not enough valid data points to calibrate. Need at least 7.');
end

% Isolate the first 7 points for linear calibration
calib_obs = obs(1:7);
calib_pred = pred_raw(1:7);

% Fit linear model: InSitu = slope * Raw_RF + intercept
p = polyfit(calib_pred, calib_obs, 1);
slope = p(1);
intercept = p(2);

fprintf('\nApplied Regional Calibration based on first 7 points:\n');
fprintf('Calibrated TSS = (%.4f * Raw_RF_TSS) + %.4f\n', slope, intercept);

% Apply calibration to all predictions
pred_calibrated = (pred_raw * slope) + intercept;
pred_calibrated = max(0.01, pred_calibrated); % Ensure positive values for log plot

%% Step 6: Statistical Analysis
RMSE = sqrt(mean((pred_calibrated - obs).^2));
R_mat = corrcoef(obs, pred_calibrated);
Rsq = R_mat(1,2)^2;

fprintf('Validation Stats (Calibrated vs In Situ):\n N = %d\n R^2 = %.4f\n RMSE = %.4f mg/L\n', ...
    N_total, Rsq, RMSE);

%% Step 7: Save Results to CSV
output_csv_name = 'Calibrated_vs_Insitu_TSS_Matchups.csv';
results_table = table((1:N_total)', obs, pred_raw, pred_calibrated, ...
    'VariableNames', {'Point_ID', 'InSitu_TSS', 'Raw_RF_TSS', 'Calibrated_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved matchup results to: %s\n', output_csv_name);

%% Step 8: Plot Calibrated vs. In Situ TSS Matchups (Log Scale)
fig = figure('Position', [100, 100, 700, 600]);

% 1. Plot all validation points (blue circles)
h_val = scatter(obs, pred_calibrated, 70, 'b', 'filled', 'MarkerEdgeColor', 'k');
hold on;

% 2. Highlight calibration subset (first 7 points in red)
h_cal = scatter(obs(1:7), pred_calibrated(1:7), 70, 'r', 'filled', 'MarkerEdgeColor', 'k');

% 3. Draw 1:1 Reference Line
min_val = min([obs; pred_calibrated]) * 0.8;
max_val = max([obs; pred_calibrated]) * 1.25;
h_line = plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1.5);

% 4. Add point numbers as labels
for i = 1:N_total
    text(obs(i) * 1.03, pred_calibrated(i), sprintf('%d', i), ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
        'FontSize', 9, 'Color', [0.2 0.2 0.2]);
end

% Set Logarithmic scale for both axes
set(gca, 'XScale', 'log', 'YScale', 'log'); 

% Formatting axes
grid on; box on;
axis([min_val max_val min_val max_val]);
xlabel('In Situ TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Calibrated Predicted TSS (mg L^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
title('PACE Matchups: Calibrated TSS vs. In Situ TSS', 'FontSize', 14);

% Add statistics box
stats_str = sprintf('N = %d\nR^2 = %.2f\nRMSE = %.2f mg L^{-1}', N_total, Rsq, RMSE);
text_x = min_val * 1.15; 
text_y = max_val * 0.70; 
text(text_x, text_y, stats_str, 'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k');

% Add legend
legend([h_val, h_cal, h_line], ...
    {'Validation Matchups', 'Calibration Subset (Points 1–7)', '1:1 Line'}, ...
    'Location', 'southeast', 'FontSize', 10);

%% Step 9: Save the Figure
output_fig_name = 'Calibrated_vs_Insitu_TSS_Matchups_LogScale.png';
try
    exportgraphics(fig, output_fig_name, 'Resolution', 300);
catch
    saveas(fig, output_fig_name);
end
fprintf('Saved scatter plot to: %s\n', output_fig_name);

%%
%% Step 1: Load Saved RF Model & Set Graphic Defaults for Q1 Journals
clc; clear; close all;

% Set publication graphic defaults (Arial font for Q1 journals)
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% File Paths
mat_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat';
csv_file = 'St_SRC_Optics_cruise_PACE_Apr2025_Oct2026.csv';

% Load TreeBagger model and wavelengths
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load Matchup CSV & Map Wavelengths
fprintf('Loading matchup data: %s\n', csv_file);
matchups = readtable(csv_file);

% Extract in-situ TSS
tss_insitu = matchups.TSS; 

% Match model wavelength columns
N_wl = length(model_wl);
expected_cols = cell(1, N_wl);
for i = 1:N_wl
    expected_cols{i} = sprintf('Rrs_%d', round(model_wl(i))); 
end

% Fallback mapping for exact wavelength match
csv_vars = matchups.Properties.VariableNames;
missing_cols = ~ismember(expected_cols, csv_vars);
if any(missing_cols)
    csv_rrs_cols = csv_vars(startsWith(csv_vars, 'Rrs_', 'IgnoreCase', true));
    csv_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), csv_rrs_cols);
    for m = find(missing_cols)
        [~, nearest_idx] = min(abs(csv_wls - model_wl(m)));
        expected_cols{m} = csv_rrs_cols{nearest_idx};
    end
end
rrs_matrix = table2array(matchups(:, expected_cols));

%% Step 3: Spectral Smoothing & 2nd Derivative
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21); % 3rd-order, 21-point Savitzky-Golay
end
first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: RF Prediction & Regional Scaling
tss_rf_raw = predict(rf_model, second_deriv_matrix);
tss_rf_raw = max(0, min(1000, tss_rf_raw));

valid_idx = ~isnan(tss_insitu) & ~isnan(tss_rf_raw);
obs = tss_insitu(valid_idx);
pred_raw = tss_rf_raw(valid_idx);
N_total = length(obs);

if N_total < 7
    error('Insufficient data points (minimum 7 required).');
end

% Regional linear fit based on first 7 points
p = polyfit(pred_raw(1:7), obs(1:7), 1);
slope = p(1); intercept = p(2);

% Final predicted TSS values
tss_predicted = max(0, (pred_raw * slope) + intercept);

%% Step 5: Statistical Calculation
RMSE = sqrt(mean((tss_predicted - obs).^2));
R_mat = corrcoef(obs, tss_predicted);
Rsq = R_mat(1,2)^2;

fprintf('\nValidation Stats:\n N = %d\n R^2 = %.4f\n RMSE = %.4f mg/L\n', ...
    N_total, Rsq, RMSE);

%% Step 6: Save Predicted vs In-Situ Matchup Data to CSV
output_csv_name = 'Predicted_vs_Insitu_TSS_Matchups_v2.csv';
results_table = table((1:N_total)', obs, tss_predicted, ...
    'VariableNames', {'Point_ID', 'InSitu_TSS', 'Predicted_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved matchup results to: %s\n', output_csv_name);

%% Step 7: Create Q1 Journal Publication Figure (Linear Scale)
fig = figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'PaperPositionMode', 'auto');
ax = axes('Parent', fig);

hold(ax, 'on');

% Colors
color_pts = [0.12 0.47 0.71];   % Deep Blue
color_line = [0.25 0.25 0.25];  % Dark Charcoal for 1:1 Line

% Calculate limits
all_vals = [obs; tss_predicted];
min_val = max(0, floor(min(all_vals) - 4));
max_val = ceil(max(all_vals) + 4);
range_span = max_val - min_val;

% 1. Draw 1:1 Line
h_line = plot(ax, [min_val, max_val], [min_val, max_val], '--', ...
    'Color', color_line, 'LineWidth', 1.2);

% 2. Scatter Points (Smaller circles: size = 35)
h_pts = scatter(ax, obs, tss_predicted, 35, color_pts, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% 3. Add Point Number Labels
for i = 1:N_total
    text(ax, obs(i) + (range_span * 0.018), tss_predicted(i), sprintf('%d', i), ...
        'FontSize', 8, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'middle');
end

% Axes & Grid Formatting
grid(ax, 'on');
box(ax, 'on');
set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 1.0, ...
    'XColor', 'k', 'YColor', 'k');

xlim(ax, [min_val max_val]);
ylim(ax, [min_val max_val]);

xlabel(ax, 'In Situ TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);
ylabel(ax, 'Predicted TSS (mg L^{-1})', 'FontWeight', 'bold', 'FontSize', 10);

% 4. Add Statistical Summary Panel (Bottom Right)
stats_str = { ...
    sprintf('N = %d', N_total), ...
    sprintf('R^2 = %.2f', Rsq), ...
    sprintf('RMSE = %.2f mg L^{-1}', RMSE) ...
};

text(ax, min_val + (range_span * 0.65), min_val + (range_span * 0.22), stats_str, ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1 0.9], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);

% 5. Legend in Top Left
legend(ax, [h_pts, h_line], ...
    {'PACE Matchups', '1:1 Line'}, ...
    'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');

hold(ax, 'off');

%% Step 8: Export Figure to Vector PDF and High-Res PNG
pdf_filename = 'Fig_PACE_TSS_Predicted_Q1_v2.pdf';
png_filename = 'Fig_PACE_TSS_Predicted_Q1_v2.png';

% Vector PDF for publication submission
exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

% High-Resolution PNG (300 DPI)
exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);



%% Going back to original model without caliration
%% Step 1: Load Saved RF Model & Set Graphic Defaults for Q1 Journals
clc; clear; close all;

% Set publication graphic defaults (Arial font for Q1 journals)
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% File Paths
% mat_file = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_30Sept25.mat';
mat_file = '/Users/masud/WorkingFile/QWIP/PACE_rf_model_2nd_deriv_best_Simul_insitu_valid_16Oct26.mat';
csv_file = 'St_cruise_PACE_26_09_2026.csv';

% Load TreeBagger model and wavelengths array
load(mat_file, 'rf_model', 'wavelengths');
model_wl = wavelengths; 

%% Step 2: Load New Matchup CSV & Extract Wavelengths
fprintf('Loading new matchup data from: %s\n', csv_file);
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
fprintf('Extracted %d Rrs columns matching the model wavelengths.\n', size(rrs_matrix, 2));

%% Step 3: Compute Smoothed Spectra and Second Derivatives
N_matchups = size(rrs_matrix, 1);
smoothed_spectra = zeros(N_matchups, length(model_wl), 'double');

fprintf('Applying Savitzky-Golay smoothing & 2nd derivative...\n');
for i = 1:N_matchups
    spec = max(rrs_matrix(i, :), eps);
    smoothed_spectra(i, :) = sgolayfilt(spec, 3, 21); % 3rd-order, 21-point Savitzky-Golay
end

first_deriv = gradient(smoothed_spectra, model_wl, 2);
second_deriv_matrix = gradient(first_deriv, model_wl, 2);

%% Step 4: Apply Original RF Model Directly (Raw Predictions)
fprintf('Predicting TSS using original RF Model...\n');
tss_predicted = predict(rf_model, second_deriv_matrix);
tss_predicted = max(0, min(1000, tss_predicted)); % Clip bounds [0, 1000]

%% Step 5: Statistical Analysis
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

fprintf('\nOriginal RF Model Evaluation Stats:\n N = %d\n R^2 = %.4f\n RMSE = %.4f mg/L\n', ...
    N_total, Rsq, RMSE);

%% Step 6: Save Results to CSV (New Output Name)
output_csv_name = 'Original_RF_vs_Insitu_TSS_Matchups_Oct2026.csv';
results_table = table((1:N_total)', obs, pred, ...
    'VariableNames', {'Point_ID', 'InSitu_TSS', 'Original_RF_Predicted_TSS'});
writetable(results_table, output_csv_name);
fprintf('Saved matchup results to: %s\n', output_csv_name);

%% Step 7: Plot Q1 Journal Publication Figure (Linear Scale)
fig = figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'PaperPositionMode', 'auto');
ax = axes('Parent', fig);

hold(ax, 'on');

% Styling parameters
color_pts = [0.12 0.47 0.71];   % Deep Blue
color_line = [0.25 0.25 0.25];  % Dark Charcoal for 1:1 Line

% Calculate plot limits
all_vals = [obs; pred];
min_val = max(0, floor(min(all_vals) - 4));
max_val = ceil(max(all_vals) + 4);
range_span = max_val - min_val;

% 1. Draw 1:1 Line
h_line = plot(ax, [min_val, max_val], [min_val, max_val], '--', ...
    'Color', color_line, 'LineWidth', 1.2);

% 2. Scatter Points (Small circles = 35)
h_pts = scatter(ax, obs, pred, 35, color_pts, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% 3. Add Point Numbers
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

% 4. Add Statistical Summary Box (Bottom Right)
stats_str = { ...
    sprintf('N = %d', N_total), ...
    sprintf('R^2 = %.2f', Rsq), ...
    sprintf('RMSE = %.2f mg L^{-1}', RMSE) ...
};

text(ax, min_val + (range_span * 0.65), min_val + (range_span * 0.22), stats_str, ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 1 0.9], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);

% 5. Add Legend in Top Left
legend(ax, [h_pts, h_line], ...
    {'PACE Matchups', '1:1 Line'}, ...
    'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');

hold(ax, 'off');

%% Step 8: Export Figure to Vector PDF and High-Res PNG (New Names)
pdf_filename = 'Fig_Original_RF_TSS_Matchups_Oct2026.pdf';
png_filename = 'Fig_Original_RF_TSS_Matchups_Oct2026.png';

% Vector PDF for publication submission
exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

% High-Resolution 300 DPI PNG
exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);
