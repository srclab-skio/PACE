%% Sorad and PACE comparison spectra for 2 cruises
%% Step 1: Initialize Workspace & Set Publication Typography Defaults
clc; clear; close all;

% Set graphic defaults to Arial for Q1 journal publication standards
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 8.5);
set(0, 'DefaultTextFontSize', 8.5);

% Define File Paths
pace_dir  = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/CleanedRrs';
sorad_dir = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';

pace_file  = fullfile(pace_dir, 'St_SRC_Optics_cruise_PACE_Apr2025_Oct2026.csv');
sorad_file = fullfile(sorad_dir, 'Optics_and_St_cruise_SoRad_Matchups_10min.csv');

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

%% Step 3: Filter Wavelengths (400 to 719 nm) & Isolate 18 Stations
pace_wl_mask  = (pace_wls >= 400) & (pace_wls <= 719);
sorad_wl_mask = (sorad_wls >= 400) & (sorad_wls <= 719);

pace_wls_sub  = pace_wls(pace_wl_mask);
sorad_wls_sub = sorad_wls(sorad_wl_mask);

N_plot = min([18, height(pace_table), height(sorad_table)]);

pace_spectra  = table2array(pace_table(1:N_plot, pace_rrs_cols));
sorad_spectra = table2array(sorad_table(1:N_plot, sorad_rrs_cols));

pace_spectra_sub  = pace_spectra(:, pace_wl_mask);
sorad_spectra_sub = sorad_spectra(:, sorad_wl_mask);

station_ids = pace_table.St(1:N_plot);

%% Step 4: Create Multi-Panel Subplot Figure (5x4 Grid for 18 Stations)
fig = figure('Units', 'centimeters', 'Position', [1, 1, 24, 26], 'PaperPositionMode', 'auto');

% Generate 18 distinct Q1 color palette lines
color_palette = lines(N_plot);

for i = 1:N_plot
    ax = subplot(5, 4, i);
    hold(ax, 'on');
    
    c = color_palette(i, :);
    
    % Plot PACE spectrum (Solid Line, width = 1.3)
    h_p = plot(ax, pace_wls_sub, pace_spectra_sub(i, :), '-', ...
        'Color', c, 'LineWidth', 1.3);
    
    % Plot SoRad spectrum (Dashed Line, width = 1.1)
    h_s = plot(ax, sorad_wls_sub, sorad_spectra_sub(i, :), '--', ...
        'Color', c * 0.7, 'LineWidth', 1.1);
    
    % Formatting Individual Axes
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'LineWidth', 0.8, ...
        'XColor', 'k', 'YColor', 'k');
    
    xlim(ax, [400 719]);
    title(ax, sprintf('Station %d', station_ids(i)), 'FontWeight', 'bold', 'FontSize', 8.5);
    
    % Add X-label conditionally to bottom row subplots
    if i >= 15
        xlabel(ax, '\lambda (nm)', 'FontWeight', 'bold', 'FontSize', 8.0);
    end
    % Add Y-label conditionally to left column subplots
    if mod(i, 4) == 1
        ylabel(ax, 'R_{rs} (sr^{-1})', 'FontWeight', 'bold', 'FontSize', 8.0);
    end
    
    % Custom Compact Legend per Subplot Panel
    legend(ax, [h_p, h_s], {'PACE', 'SoRad'}, ...
        'Location', 'northeast', 'FontSize', 6.5, 'Box', 'on');
    
    hold(ax, 'off');
end

% Hide unused 19th and 20th subplot panels
ax_unused19 = subplot(5, 4, 19); axis(ax_unused19, 'off');
ax_unused20 = subplot(5, 4, 20); axis(ax_unused20, 'off');

%% Step 5: Export Vector PDF & High-Res PNG with Identifiable Filenames
pdf_filename = fullfile(sorad_dir, 'Fig_PACE_vs_SoRad_Spectra_18Stations_Subplots_Q1.pdf');
png_filename = fullfile(sorad_dir, 'Fig_PACE_vs_SoRad_Spectra_18Stations_Subplots_Q1.png');

% Export vector PDF for Q1 manuscript submission
exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
fprintf('Exported Vector PDF: %s\n', pdf_filename);

% Export 300 DPI PNG
exportgraphics(fig, png_filename, 'Resolution', 300);
fprintf('Exported 300 DPI PNG: %s\n', png_filename);

