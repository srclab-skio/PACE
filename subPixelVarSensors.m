%% GOK
% Clear workspace and command window
clear;
clc;

% Load the .mat file
matFilePath = '/Users/masud/OneDriveUGA/CruiseRVSavannah/Subpixel/updated_Rrs_SoRad_PACE_OLCI.mat';
data = load(matFilePath);
RrsStruct = data.Rrs_SoRad_PACE_OLCI;
% Define stations
stations = {'St10', 'St11', 'St12', 'St14', 'St15', 'St16'};

% Initialize metrics storage
metrics = struct('Station', {}, 'SoRad_vs_PACE_RMSE', {}, 'SoRad_vs_OLCI_RMSE', {}, ...
                 'PACE_vs_OLCI_RMSE', {}, 'SoRad_vs_PACE_MAPD', {}, 'SoRad_vs_OLCI_MAPD', {}, ...
                 'PACE_vs_OLCI_MAPD', {}, 'SoRad_vs_PACE_SA', {}, 'SoRad_vs_OLCI_SA', {}, ...
                 'PACE_vs_OLCI_SA', {});

%% For plots only
% Comparison plots
% Comparison plots
figure('Name', 'Rrs Comparison by Station', 'Position', [100, 100, 1200, 800]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % No space between subplots

for i = 1:length(stations)
    stName = stations{i};
    if ~isfield(RrsStruct.Stations, stName)
        warning('Station %s not found in Stations structure', stName);
        continue;
    end
    stData = RrsStruct.Stations.(stName);
    
    % Compute mean and std for each sensor
    soradMean = mean(stData.Sorad, 1, 'omitnan');
    soradStd = std(stData.Sorad, 0, 1, 'omitnan');
    paceMean = mean(stData.PACE, 1, 'omitnan');
    paceStd = std(stData.PACE, 0, 1, 'omitnan');
    olciMean = mean(stData.OLCI, 1, 'omitnan');
    olciStd = std(stData.OLCI, 0, 1, 'omitnan');
    
    % Subplot for each station
    nexttile;
    hold on;
    
    % Plot SoRad with error bars (if multiple rows)
    if size(stData.Sorad, 1) > 1
        errorbar(RrsStruct.wl_Sorad, soradMean, soradStd, 'b-', 'DisplayName', 'SoRad','LineWidth', 1);
    else
        plot(RrsStruct.wl_Sorad, soradMean, 'b-', 'DisplayName', 'SoRad','LineWidth', 1);
    end
    
    % Interpolate PACE to SoRad wavelengths (400-800 nm)
    idx_wl = RrsStruct.wl_Sorad >= 400 & RrsStruct.wl_Sorad <= 800;
    wl_common = RrsStruct.wl_Sorad(idx_wl);
    paceMean_interp = interp1(RrsStruct.wl_PACE, paceMean, wl_common, 'linear', 'extrap');
    paceStd_interp = interp1(RrsStruct.wl_PACE, paceStd, wl_common, 'linear', 'extrap');
    if size(stData.PACE, 1) > 1
        errorbar(RrsStruct.wl_PACE, paceMean, paceStd, 'r-', 'DisplayName', 'PACE','LineWidth', 1);
    else
        plot(RrsStruct.wl_PACE, paceMean, 'r-', 'DisplayName', 'PACE','LineWidth', 1);
    end
    
    % Interpolate OLCI to SoRad wavelengths
    olciMean_interp = interp1(RrsStruct.wl_OLCI, olciMean, wl_common, 'linear', 'extrap');
    olciStd_interp = interp1(RrsStruct.wl_OLCI, olciStd, wl_common, 'linear', 'extrap');
    if size(stData.OLCI, 1) > 1
        errorbar(RrsStruct.wl_OLCI, olciMean, olciStd, 'g-', 'DisplayName', 'OLCI','LineWidth', 1);
    else
        plot(RrsStruct.wl_OLCI, olciMean, 'g-', 'DisplayName', 'OLCI','LineWidth', 1);
    end
    
    % Plot settings
    xlabel('Wavelength (nm)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
    ylabel('Rrs (sr^{-1})', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
    title(sprintf('Station %s', stName), 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial');
    
    % Display number of spectra for OLCI and SoRad on the figure
    text(0.5, 0.95, sprintf('SoRad Spectra: %d, OLCI Spectra: %d', size(stData.Sorad, 1), size(stData.OLCI, 1)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 12);
    
    % Add the station number inside the figure area
    % text(0.05, 0.9, sprintf('Station: %d', i), 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial');
    
    % Set grid and axis properties
    legend('show', 'FontWeight', 'bold', 'FontSize', 10,'Location','best');
    grid off;
    xlim([400 800]);
    ylim([-0.001 0.025]);
    
    % Add outline for each subplot
    set(gca, 'Box', 'on', 'LineWidth', 2);  % Make the box thicker
    
    hold off;
end

% Save the figure as PNG
saveas(gcf, 'Rrs_Comparison_by_Station.png');

filename = fullfile('/Users/masud/OneDriveUGA/CruiseRVSavannah/Comparison_plots/', "Good_6_with_Error_bars" + ""  + ".png");
exportgraphics(gcf,filename,'Resolution',600)







%%

%% GOK
%% GOK
% Clear workspace and command window
clear;
clc;

% Load the .mat file
matFilePath = '/Users/masud/OneDriveUGA/CruiseRVSavannah/Subpixel/updated_Rrs_SoRad_PACE_OLCI.mat';
data = load(matFilePath);
RrsStruct = data.Rrs_SoRad_PACE_OLCI;

% Define stations
stations = {'St10', 'St11', 'St12', 'St14', 'St15', 'St16'};

% Initialize metrics storage
metrics = struct('Station', {}, 'SoRad_vs_PACE_RMSE', {}, 'SoRad_vs_OLCI_RMSE', {}, ...
                 'PACE_vs_OLCI_RMSE', {}, 'SoRad_vs_PACE_MAPD', {}, 'SoRad_vs_OLCI_MAPD', {}, ...
                 'PACE_vs_OLCI_MAPD', {}, 'SoRad_vs_PACE_SA', {}, 'SoRad_vs_OLCI_SA', {}, ...
                 'PACE_vs_OLCI_SA', {});

% Comparison plots
figure('Name', 'Rrs Comparison by Station', 'Position', [100, 100, 1200, 800]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(stations)
    stName = stations{i};
    if ~isfield(RrsStruct.Stations, stName)
        warning('Station %s not found in Stations structure', stName);
        continue;
    end
    stData = RrsStruct.Stations.(stName);
    
    % Compute mean and std for each sensor
    soradMean = mean(stData.Sorad, 1, 'omitnan');
    soradStd = std(stData.Sorad, 0, 1, 'omitnan');
    paceMean = mean(stData.PACE, 1, 'omitnan');
    paceStd = std(stData.PACE, 0, 1, 'omitnan');
    olciMean = mean(stData.OLCI, 1, 'omitnan');
    olciStd = std(stData.OLCI, 0, 1, 'omitnan');
    
    % Subplot for each station
    nexttile;
    hold on;
    
    % Plot SoRad with error bars (if multiple rows) - thicker lines
    if size(stData.Sorad, 1) > 1
        errorbar(RrsStruct.wl_Sorad, soradMean, soradStd, 'b-', 'LineWidth', 2, 'DisplayName', 'SoRad');
    else
        plot(RrsStruct.wl_Sorad, soradMean, 'b-', 'LineWidth', 2, 'DisplayName', 'SoRad');
    end
    
    % Interpolate PACE to SoRad wavelengths (400-800 nm)
    idx_wl = RrsStruct.wl_Sorad >= 400 & RrsStruct.wl_Sorad <= 800;
    wl_common = RrsStruct.wl_Sorad(idx_wl);
    paceMean_interp = interp1(RrsStruct.wl_PACE, paceMean, wl_common, 'linear', 'extrap');
    paceStd_interp = interp1(RrsStruct.wl_PACE, paceStd, wl_common, 'linear', 'extrap');
    if size(stData.PACE, 1) > 1
        errorbar(wl_common, paceMean_interp, paceStd_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'PACE');
    else
        plot(wl_common, paceMean_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'PACE');
    end
    
    % Interpolate OLCI to SoRad wavelengths
    olciMean_interp = interp1(RrsStruct.wl_OLCI, olciMean, wl_common, 'linear', 'extrap');
    olciStd_interp = interp1(RrsStruct.wl_OLCI, olciStd, wl_common, 'linear', 'extrap');
    if size(stData.OLCI, 1) > 1
        errorbar(wl_common, olciMean_interp, olciStd_interp, 'g-', 'LineWidth', 2, 'DisplayName', 'OLCI');
    else
        plot(wl_common, olciMean_interp, 'g-', 'LineWidth', 2, 'DisplayName', 'OLCI');
    end
    
    % Plot settings with bold fonts
    xlabel('Wavelength (nm)', 'FontWeight', 'bold');
    ylabel('Rrs (sr^{-1})', 'FontWeight', 'bold');
    title(sprintf('Station %s', stName), 'FontWeight', 'bold');
    legend('show', 'FontWeight', 'bold');
    grid on;
    xlim([400 800]);
    hold off;
    
    % Calculate metrics (SoRad as reference wavelength grid)
    metrics(i).Station = stName;
    
    % RMSE
    metrics(i).SoRad_vs_PACE_RMSE = sqrt(mean((soradMean(idx_wl) - paceMean_interp).^2, 'omitnan'));
    metrics(i).SoRad_vs_OLCI_RMSE = sqrt(mean((soradMean(idx_wl) - olciMean_interp).^2, 'omitnan'));
    metrics(i).PACE_vs_OLCI_RMSE = sqrt(mean((paceMean_interp - olciMean_interp).^2, 'omitnan'));
    
    % MAPD (avoid division by zero)
    soradMean_valid = soradMean(idx_wl);
    paceMean_valid = paceMean_interp;
    olciMean_valid = olciMean_interp;
    valid_idx = abs(soradMean_valid) > 1e-6; % Avoid division by near-zero
    metrics(i).SoRad_vs_PACE_MAPD = mean(abs((soradMean_valid(valid_idx) - paceMean_valid(valid_idx)) ./ ...
        soradMean_valid(valid_idx)) * 100, 'omitnan');
    valid_idx = abs(soradMean_valid) > 1e-6;
    metrics(i).SoRad_vs_OLCI_MAPD = mean(abs((soradMean_valid(valid_idx) - olciMean_valid(valid_idx)) ./ ...
        soradMean_valid(valid_idx)) * 100, 'omitnan');
    valid_idx = abs(paceMean_valid) > 1e-6;
    metrics(i).PACE_vs_OLCI_MAPD = mean(abs((paceMean_valid(valid_idx) - olciMean_valid(valid_idx)) ./ ...
        paceMean_valid(valid_idx)) * 100, 'omitnan');
    
    % Spectral Angle (in degrees)
    dot_product = @(a, b) sum(a .* b, 'omitnan');
    norm_vector = @(a) sqrt(sum(a.^2, 'omitnan'));
    sa = @(a, b) acosd(min(dot_product(a, b) / (norm_vector(a) * norm_vector(b)), 1));
    metrics(i).SoRad_vs_PACE_SA = sa(soradMean(idx_wl), paceMean_interp);
    metrics(i).SoRad_vs_OLCI_SA = sa(soradMean(idx_wl), olciMean_interp);
    metrics(i).PACE_vs_OLCI_SA = sa(paceMean_interp, olciMean_interp);
end

% Display metrics
disp('Metrics for each station:');
for i = 1:length(metrics)
    fprintf('Station %s:\n', metrics(i).Station);
    fprintf('  SoRad vs PACE RMSE: %.6f sr^-1\n', metrics(i).SoRad_vs_PACE_RMSE);
    fprintf('  SoRad vs OLCI RMSE: %.6f sr^-1\n', metrics(i).SoRad_vs_OLCI_RMSE);
    fprintf('  PACE vs OLCI RMSE: %.6f sr^-1\n', metrics(i).PACE_vs_OLCI_RMSE);
    fprintf('  SoRad vs PACE MAPD: %.2f%%\n', metrics(i).SoRad_vs_PACE_MAPD);
    fprintf('  SoRad vs OLCI MAPD: %.2f%%\n', metrics(i).SoRad_vs_OLCI_MAPD);
    fprintf('  PACE vs OLCI MAPD: %.2f%%\n', metrics(i).PACE_vs_OLCI_MAPD);
    fprintf('  SoRad vs PACE SA: %.2f degrees\n', metrics(i).SoRad_vs_PACE_SA);
    fprintf('  SoRad vs OLCI SA: %.2f degrees\n', metrics(i).SoRad_vs_OLCI_SA);
    fprintf('  PACE vs OLCI SA: %.2f degrees\n', metrics(i).PACE_vs_OLCI_SA);
end

% Export metrics to separate CSV files
% Get current date and time
currentTime = datestr(now, 'yyyy-mm-dd_HHMM');
baseDir = fileparts(matFilePath);

% RMSE CSV
rmseFilePath = fullfile(baseDir, sprintf('rmse_output_%s.csv', currentTime));
rmseData = struct2table(metrics);
rmseFields = {'Station', 'SoRad_vs_PACE_RMSE', 'SoRad_vs_OLCI_RMSE', 'PACE_vs_OLCI_RMSE'};
rmseData = rmseData(:, rmseFields);
writetable(rmseData, rmseFilePath);
disp(['RMSE metrics saved to ' rmseFilePath]);

% MAPD CSV
mapdFilePath = fullfile(baseDir, sprintf('mapd_output_%s.csv', currentTime));
mapdData = struct2table(metrics);
mapdFields = {'Station', 'SoRad_vs_PACE_MAPD', 'SoRad_vs_OLCI_MAPD', 'PACE_vs_OLCI_MAPD'};
mapdData = mapdData(:, mapdFields);
writetable(mapdData, mapdFilePath);
disp(['MAPD metrics saved to ' mapdFilePath]);

% Spectral Angle CSV
saFilePath = fullfile(baseDir, sprintf('sa_output_%s.csv', currentTime));
saData = struct2table(metrics);
saFields = {'Station', 'SoRad_vs_PACE_SA', 'SoRad_vs_OLCI_SA', 'PACE_vs_OLCI_SA'};
saData = saData(:, saFields);
writetable(saData, saFilePath);
disp(['Spectral Angle metrics saved to ' saFilePath]);

%%




%% GOK
% Clear workspace and command window
clear;
clc;

% Load the .mat file
matFilePath = '/Users/masud/OneDriveUGA/CruiseRVSavannah/Subpixel/updated_Rrs_SoRad_PACE_OLCI.mat';
data = load(matFilePath);
RrsStruct = data.Rrs_SoRad_PACE_OLCI;

% Define stations
stations = {'St10', 'St11', 'St12', 'St14', 'St15', 'St16'};

% Initialize metrics storage
metrics = struct('Station', {}, 'SoRad_vs_PACE_RMSE', {}, 'SoRad_vs_OLCI_RMSE', {}, ...
                 'PACE_vs_OLCI_RMSE', {}, 'SoRad_vs_PACE_MAPD', {}, 'SoRad_vs_OLCI_MAPD', {}, ...
                 'PACE_vs_OLCI_MAPD', {}, 'SoRad_vs_PACE_SA', {}, 'SoRad_vs_OLCI_SA', {}, ...
                 'PACE_vs_OLCI_SA', {});

% Comparison plots
figure('Name', 'Rrs Comparison by Station', 'Position', [100, 100, 1200, 800]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(stations)
    stName = stations{i};
    if ~isfield(RrsStruct.Stations, stName)
        warning('Station %s not found in Stations structure', stName);
        continue;
    end
    stData = RrsStruct.Stations.(stName);
    
    % Compute mean and std for each sensor
    soradMean = mean(stData.Sorad, 1, 'omitnan');
    soradStd = std(stData.Sorad, 0, 1, 'omitnan');
    paceMean = mean(stData.PACE, 1, 'omitnan');
    paceStd = std(stData.PACE, 0, 1, 'omitnan');
    olciMean = mean(stData.OLCI, 1, 'omitnan');
    olciStd = std(stData.OLCI, 0, 1, 'omitnan');
    
    % Subplot for each station
    nexttile;
    hold on;
    
    % Plot SoRad with error bars (if multiple rows) - thicker lines
    if size(stData.Sorad, 1) > 1
        errorbar(RrsStruct.wl_Sorad, soradMean, soradStd, 'b-', 'LineWidth', 2, 'DisplayName', 'SoRad');
    else
        plot(RrsStruct.wl_Sorad, soradMean, 'b-', 'LineWidth', 2, 'DisplayName', 'SoRad');
    end
    
    % Interpolate PACE to SoRad wavelengths (400-800 nm)
    idx_wl = RrsStruct.wl_Sorad >= 400 & RrsStruct.wl_Sorad <= 800;
    wl_common = RrsStruct.wl_Sorad(idx_wl);
    paceMean_interp = interp1(RrsStruct.wl_PACE, paceMean, wl_common, 'linear', 'extrap');
    paceStd_interp = interp1(RrsStruct.wl_PACE, paceStd, wl_common, 'linear', 'extrap');
    if size(stData.PACE, 1) > 1
        errorbar(wl_common, paceMean_interp, paceStd_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'PACE');
    else
        plot(wl_common, paceMean_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'PACE');
    end
    
    % Interpolate OLCI to SoRad wavelengths
    olciMean_interp = interp1(RrsStruct.wl_OLCI, olciMean, wl_common, 'linear', 'extrap');
    olciStd_interp = interp1(RrsStruct.wl_OLCI, olciStd, wl_common, 'linear', 'extrap');
    if size(stData.OLCI, 1) > 1
        errorbar(wl_common, olciMean_interp, olciStd_interp, 'g-', 'LineWidth', 2, 'DisplayName', 'OLCI');
    else
        plot(wl_common, olciMean_interp, 'g-', 'LineWidth', 2, 'DisplayName', 'OLCI');
    end
    
    % Plot settings with bold fonts
    xlabel('Wavelength (nm)', 'FontWeight', 'bold');
    ylabel('Rrs (sr^{-1})', 'FontWeight', 'bold');
    title(sprintf('Station %s', stName), 'FontWeight', 'bold');
    legend('show', 'FontWeight', 'bold');
    grid on;
    xlim([400 800]);
    hold off;
    
    % Calculate metrics (SoRad as reference wavelength grid)
    metrics(i).Station = stName;
    
    % RMSE
    metrics(i).SoRad_vs_PACE_RMSE = sqrt(mean((soradMean(idx_wl) - paceMean_interp).^2, 'omitnan'));
    metrics(i).SoRad_vs_OLCI_RMSE = sqrt(mean((soradMean(idx_wl) - olciMean_interp).^2, 'omitnan'));
    metrics(i).PACE_vs_OLCI_RMSE = sqrt(mean((paceMean_interp - olciMean_interp).^2, 'omitnan'));
    
    % MAPD (avoid division by zero)
    soradMean_valid = soradMean(idx_wl);
    paceMean_valid = paceMean_interp;
    olciMean_valid = olciMean_interp;
    valid_idx = abs(soradMean_valid) > 1e-6; % Avoid division by near-zero
    metrics(i).SoRad_vs_PACE_MAPD = mean(abs((soradMean_valid(valid_idx) - paceMean_valid(valid_idx)) ./ ...
        soradMean_valid(valid_idx)) * 100, 'omitnan');
    valid_idx = abs(soradMean_valid) > 1e-6;
    metrics(i).SoRad_vs_OLCI_MAPD = mean(abs((soradMean_valid(valid_idx) - olciMean_valid(valid_idx)) ./ ...
        soradMean_valid(valid_idx)) * 100, 'omitnan');
    valid_idx = abs(paceMean_valid) > 1e-6;
    metrics(i).PACE_vs_OLCI_MAPD = mean(abs((paceMean_valid(valid_idx) - olciMean_valid(valid_idx)) ./ ...
        paceMean_valid(valid_idx)) * 100, 'omitnan');
    
    % Spectral Angle (in degrees)
    dot_product = @(a, b) sum(a .* b, 'omitnan');
    norm_vector = @(a) sqrt(sum(a.^2, 'omitnan'));
    sa = @(a, b) acosd(min(dot_product(a, b) / (norm_vector(a) * norm_vector(b)), 1));
    metrics(i).SoRad_vs_PACE_SA = sa(soradMean(idx_wl), paceMean_interp);
    metrics(i).SoRad_vs_OLCI_SA = sa(soradMean(idx_wl), olciMean_interp);
    metrics(i).PACE_vs_OLCI_SA = sa(paceMean_interp, olciMean_interp);
end

% Display metrics
disp('Metrics for each station:');
for i = 1:length(metrics)
    fprintf('Station %s:\n', metrics(i).Station);
    fprintf('  SoRad vs PACE RMSE: %.6f sr^-1\n', metrics(i).SoRad_vs_PACE_RMSE);
    fprintf('  SoRad vs OLCI RMSE: %.6f sr^-1\n', metrics(i).SoRad_vs_OLCI_RMSE);
    fprintf('  PACE vs OLCI RMSE: %.6f sr^-1\n', metrics(i).PACE_vs_OLCI_RMSE);
    fprintf('  SoRad vs PACE MAPD: %.2f%%\n', metrics(i).SoRad_vs_PACE_MAPD);
    fprintf('  SoRad vs OLCI MAPD: %.2f%%\n', metrics(i).SoRad_vs_OLCI_MAPD);
    fprintf('  PACE vs OLCI MAPD: %.2f%%\n', metrics(i).PACE_vs_OLCI_MAPD);
    fprintf('  SoRad vs PACE SA: %.2f degrees\n', metrics(i).SoRad_vs_PACE_SA);
    fprintf('  SoRad vs OLCI SA: %.2f degrees\n', metrics(i).SoRad_vs_OLCI_SA);
    fprintf('  PACE vs OLCI SA: %.2f degrees\n', metrics(i).PACE_vs_OLCI_SA);
end

% Export metrics to .mat file
% Get current date and time
currentTime = datestr(now, 'yyyy-mm-dd_HHMM');
matFilePathOut = fullfile(fileparts(matFilePath), sprintf('metrics_output_%s.mat', currentTime));

% Save metrics structure to .mat file
save(matFilePathOut, 'metrics');
disp(['Metrics saved to ' matFilePathOut]);
%%
%% GOK
% Clear workspace and command window
clear;
clc;

% Load the .mat file
matFilePath = '/Users/masud/OneDriveUGA/CruiseRVSavannah/Subpixel/updated_Rrs_SoRad_PACE_OLCI.mat';
data = load(matFilePath);
RrsStruct = data.Rrs_SoRad_PACE_OLCI;

% Define stations
stations = {'St10', 'St11', 'St12', 'St14', 'St15', 'St16'};

% Initialize metrics storage
metrics = struct('Station', {}, 'SoRad_vs_PACE_RMSE', {}, 'SoRad_vs_OLCI_RMSE', {}, ...
                 'PACE_vs_OLCI_RMSE', {}, 'SoRad_vs_PACE_MAPD', {}, 'SoRad_vs_OLCI_MAPD', {}, ...
                 'PACE_vs_OLCI_MAPD', {}, 'SoRad_vs_PACE_SA', {}, 'SoRad_vs_OLCI_SA', {}, ...
                 'PACE_vs_OLCI_SA', {});

% Calculate metrics
for i = 1:length(stations)
    stName = stations{i};
    if ~isfield(RrsStruct.Stations, stName)
        warning('Station %s not found in Stations structure', stName);
        continue;
    end
    stData = RrsStruct.Stations.(stName);
    
    % Compute mean spectra
    soradMean = mean(stData.Sorad, 1, 'omitnan');
    paceMean = mean(stData.PACE, 1, 'omitnan');
    olciMean = mean(stData.OLCI, 1, 'omitnan');
    
    % Interpolate to common wavelength grid (400-800 nm)
    idx_wl = RrsStruct.wl_Sorad >= 400 & RrsStruct.wl_Sorad <= 800;
    wl_common = RrsStruct.wl_Sorad(idx_wl);
    paceMean_interp = interp1(RrsStruct.wl_PACE, paceMean, wl_common, 'linear', 'extrap');
    olciMean_interp = interp1(RrsStruct.wl_OLCI, olciMean, wl_common, 'linear', 'extrap');
    
    % Store station name
    metrics(i).Station = stName;
    
    % RMSE
    metrics(i).SoRad_vs_PACE_RMSE = sqrt(mean((soradMean(idx_wl) - paceMean_interp).^2, 'omitnan'));
    metrics(i).SoRad_vs_OLCI_RMSE = sqrt(mean((soradMean(idx_wl) - olciMean_interp).^2, 'omitnan'));
    metrics(i).PACE_vs_OLCI_RMSE = sqrt(mean((paceMean_interp - olciMean_interp).^2, 'omitnan'));
    
    % MAPD (avoid division by zero)
    soradMean_valid = soradMean(idx_wl);
    paceMean_valid = paceMean_interp;
    olciMean_valid = olciMean_interp;
    valid_idx = abs(soradMean_valid) > 1e-6;
    metrics(i).SoRad_vs_PACE_MAPD = mean(abs((soradMean_valid(valid_idx) - paceMean_valid(valid_idx)) ./ ...
        soradMean_valid(valid_idx)) * 100, 'omitnan');
    valid_idx = abs(soradMean_valid) > 1e-6;
    metrics(i).SoRad_vs_OLCI_MAPD = mean(abs((soradMean_valid(valid_idx) - olciMean_valid(valid_idx)) ./ ...
        soradMean_valid(valid_idx)) * 100, 'omitnan');
    valid_idx = abs(paceMean_valid) > 1e-6;
    metrics(i).PACE_vs_OLCI_MAPD = mean(abs((paceMean_valid(valid_idx) - olciMean_valid(valid_idx)) ./ ...
        paceMean_valid(valid_idx)) * 100, 'omitnan');
    
    % Spectral Angle (in degrees) - not plotted but calculated for .mat
    dot_product = @(a, b) sum(a .* b, 'omitnan');
    norm_vector = @(a) sqrt(sum(a.^2, 'omitnan'));
    sa = @(a, b) acosd(min(dot_product(a, b) / (norm_vector(a) * norm_vector(b)), 1));
    metrics(i).SoRad_vs_PACE_SA = sa(soradMean(idx_wl), paceMean_interp);
    metrics(i).SoRad_vs_OLCI_SA = sa(soradMean(idx_wl), olciMean_interp);
    metrics(i).PACE_vs_OLCI_SA = sa(paceMean_interp, olciMean_interp);
end

% Display metrics
disp('Metrics for each station:');
for i = 1:length(metrics)
    fprintf('Station %s:\n', metrics(i).Station);
    fprintf('  SoRad vs PACE RMSE: %.6f sr^-1\n', metrics(i).SoRad_vs_PACE_RMSE);
    fprintf('  SoRad vs OLCI RMSE: %.6f sr^-1\n', metrics(i).SoRad_vs_OLCI_RMSE);
    fprintf('  PACE vs OLCI RMSE: %.6f sr^-1\n', metrics(i).PACE_vs_OLCI_RMSE);
    fprintf('  SoRad vs PACE MAPD: %.2f%%\n', metrics(i).SoRad_vs_PACE_MAPD);
    fprintf('  SoRad vs OLCI MAPD: %.2f%%\n', metrics(i).SoRad_vs_OLCI_MAPD);
    fprintf('  PACE vs OLCI MAPD: %.2f%%\n', metrics(i).PACE_vs_OLCI_MAPD);
    fprintf('  SoRad vs PACE SA: %.2f degrees\n', metrics(i).SoRad_vs_PACE_SA);
    fprintf('  SoRad vs OLCI SA: %.2f degrees\n', metrics(i).SoRad_vs_OLCI_SA);
    fprintf('  PACE vs OLCI SA: %.2f degrees\n', metrics(i).PACE_vs_OLCI_SA);
end

% Plot RMSE and MAPD in 6x2 subplot
figure('Name', 'RMSE and MAPD by Station', 'Position', [100, 100, 1200, 800]);
tiledlayout(6, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Define colors for sensors and markers for stations
colors = {'b', 'r', 'g'}; % Blue for SoRad vs. PACE, Red for SoRad vs. OLCI, Green for PACE vs. OLCI
markers = {'o', 's', 'd', '^', 'v', '>'}; % Circle, Square, Diamond, Triangle Up, Triangle Down, Triangle Right

for i = 1:length(stations)
    if ~isfield(metrics(i), 'Station') || isempty(metrics(i).Station)
        continue;
    end
    
    % RMSE subplot
    nexttile(i); % First column for RMSE
    hold on;
    plot(1, metrics(i).SoRad_vs_PACE_RMSE, [colors{1} markers{i}], 'LineWidth', 2, 'DisplayName', 'SoRad vs PACE');
    plot(2, metrics(i).SoRad_vs_OLCI_RMSE, [colors{2} markers{i}], 'LineWidth', 2, 'DisplayName', 'SoRad vs OLCI');
    plot(3, metrics(i).PACE_vs_OLCI_RMSE, [colors{3} markers{i}], 'LineWidth', 2, 'DisplayName', 'PACE vs OLCI');
    xlabel('Pair', 'FontWeight', 'bold');
    ylabel('RMSE (sr^{-1})', 'FontWeight', 'bold');
    title(sprintf('Station %s RMSE', metrics(i).Station), 'FontWeight', 'bold');
    legend('show', 'FontWeight', 'bold', 'Location', 'best');
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels({'SoRad vs PACE', 'SoRad vs OLCI', 'PACE vs OLCI'});
    grid on;
    hold off;
    
    % MAPD subplot
    nexttile(i + length(stations)); % Second column for MAPD
    hold on;
    plot(1, metrics(i).SoRad_vs_PACE_MAPD, [colors{1} markers{i}], 'LineWidth', 2, 'DisplayName', 'SoRad vs PACE');
    plot(2, metrics(i).SoRad_vs_OLCI_MAPD, [colors{2} markers{i}], 'LineWidth', 2, 'DisplayName', 'SoRad vs OLCI');
    plot(3, metrics(i).PACE_vs_OLCI_MAPD, [colors{3} markers{i}], 'LineWidth', 2, 'DisplayName', 'PACE vs OLCI');
    xlabel('Pair', 'FontWeight', 'bold');
    ylabel('MAPD (%)', 'FontWeight', 'bold');
    title(sprintf('Station %s MAPD', metrics(i).Station), 'FontWeight', 'bold');
    legend('show', 'FontWeight', 'bold', 'Location', 'best');
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels({'SoRad vs PACE', 'SoRad vs OLCI', 'PACE vs OLCI'});
    grid on;
    hold off;
end

% Export metrics to .mat file
% Get current date and time (06:53 PM EDT, July 15, 2025)
currentTime = datestr(now, 'yyyy-mm-dd_HHMM');
matFilePathOut = fullfile(fileparts(matFilePath), sprintf('metrics_output_%s.mat', currentTime));

% Save metrics structure to .mat file
save(matFilePathOut, 'metrics');
disp(['Metrics saved to ' matFilePathOut]);

%% using chat gpt
% Load the .mat file
load('updated_Rrs_SoRad_PACE_OLCI.mat');

% Get the station names
stations = fieldnames(Rrs_SoRad_PACE_OLCI.Stations);

% Prepare for plots
num_stations = numel(stations);

% Initialize the figure for comparison plots (2 rows, 3 columns)
figure;

% Common wavelength grid (based on the range of all wavelengths)
wl_common = 400:1:800; % Adjust as needed

% Initialize arrays for spectral angles
spectral_angles = zeros(num_stations, 3);

% Loop through each station
for i = 1:num_stations
    station = stations{i};
    
    % Extract Rrs data and wavelengths for the station
    Rrs_Sorad = Rrs_SoRad_PACE_OLCI.Stations.(station).Sorad; % [13x560]
    Rrs_PACE = Rrs_SoRad_PACE_OLCI.Stations.(station).PACE; % [1x14]
    Rrs_OLCI = Rrs_SoRad_PACE_OLCI.Stations.(station).OLCI; % [9x16]
    
    % Extract wavelengths for plotting
    wl_PACE = Rrs_SoRad_PACE_OLCI.wl_PACE; % [172x1]
    wl_Sorad = Rrs_SoRad_PACE_OLCI.wl_Sorad; % [560x1]
    wl_OLCI = Rrs_SoRad_PACE_OLCI.wl_OLCI; % [16x1]
    
    % Interpolate each Rrs to the common wavelength grid
    Rrs_Sorad_interp = interp1(wl_Sorad, Rrs_Sorad(1, :), wl_common, 'linear', 'extrap');
    Rrs_PACE_interp = interp1(wl_PACE, Rrs_PACE, wl_common, 'linear', 'extrap');
    Rrs_OLCI_interp = interp1(wl_OLCI, mean(Rrs_OLCI, 1), wl_common, 'linear', 'extrap');
    
    % Plot comparison: SoRad vs PACE vs OLCI
    subplot(2, 3, i); 
    hold on;
%     plot(wl_common, Rrs_Sorad_interp, 'r', 'LineWidth', 2); % SoRad
%     plot(wl_common, Rrs_PACE_interp, 'b', 'LineWidth', 2); % PACE
%     plot(wl_common, Rrs_OLCI_interp, 'g', 'LineWidth', 2); % OLCI
    plot(wl_common, Rrs_Sorad_interp, 'r', 'LineWidth', 2); % SoRad
    plot(wl_common, Rrs_PACE_interp, 'b', 'LineWidth', 2); % PACE
    plot(wl_common, Rrs_OLCI_interp, 'g', 'LineWidth', 2); % OLCI
    title(['', station], 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Rrs (Sr^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([400 720]);
    ylim([0.000 0.02]   )
    legend('SoRad', 'PACE', 'OLCI', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold'); % Bold and thick axis labels
    hold off;
    
    % Calculate RMSE, MAPD, and Spectral Angle for the comparison
    RMSE_SoRad_PACE = sqrt(mean((Rrs_Sorad_interp - Rrs_PACE_interp).^2)); % RMSE between SoRad and PACE
    RMSE_PACE_OLCI = sqrt(mean((Rrs_PACE_interp - Rrs_OLCI_interp).^2)); % RMSE between PACE and OLCI
    MAPD_SoRad_PACE = mean(abs(Rrs_Sorad_interp - Rrs_PACE_interp) ./ Rrs_Sorad_interp) * 100; % MAPD between SoRad and PACE
    MAPD_PACE_OLCI = mean(abs(Rrs_PACE_interp - Rrs_OLCI_interp) ./ Rrs_PACE_interp) * 100; % MAPD between PACE and OLCI
    
    % Calculate Spectral Angle between SoRad and PACE
    spectral_angle_SoRad_PACE = acos(sum(Rrs_Sorad_interp .* Rrs_PACE_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_PACE_interp))) * (180 / pi); 
    
    % Spectral Angle between PACE and OLCI
    spectral_angle_PACE_OLCI = acos(sum(Rrs_PACE_interp .* Rrs_OLCI_interp) / (norm(Rrs_PACE_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);

    % Spectral Angle between SoRad and OLCI
    spectral_angle_SoRad_OLCI = acos(sum(Rrs_Sorad_interp .* Rrs_OLCI_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);
    
    % Store Spectral Angle results for the bar plot
    spectral_angles(i, :) = [spectral_angle_SoRad_PACE, spectral_angle_PACE_OLCI, spectral_angle_SoRad_OLCI];
end

% Create a bar plot for the spectral angles
figure;
bar(spectral_angles);
set(gca, 'XTickLabel', stations, 'FontSize', 12, 'FontWeight', 'bold');
title('Spectral Angle Comparison Between Sensors', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Stations', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spectral Angle (Degrees)', 'FontSize', 12, 'FontWeight', 'bold');
legend('SoRad vs PACE', 'PACE vs OLCI', 'SoRad vs OLCI', 'FontSize', 12, 'FontWeight', 'bold');

% Save the spectral angle matrices as a CSV file
writetable(array2table(spectral_angles, 'VariableNames', {'SoRad_vs_PACE', 'PACE_vs_OLCI', 'SoRad_vs_OLCI'}), 'spectral_angles.csv');

%%
% Load the .mat file
load('updated_Rrs_SoRad_PACE_OLCI.mat');

% Get the station names
stations = fieldnames(Rrs_SoRad_PACE_OLCI.Stations);

% Prepare for plots
num_stations = numel(stations);

% Initialize arrays for spectral angles, RMSE, and MAPD
spectral_angles = zeros(num_stations, 3);
RMSE_values = zeros(num_stations, 2); % SoRad vs PACE and PACE vs OLCI
MAPD_values = zeros(num_stations, 2); % SoRad vs PACE and PACE vs OLCI

% Loop through each station
for i = 1:num_stations
    station = stations{i};
    
    % Extract Rrs data and wavelengths for the station
    Rrs_Sorad = Rrs_SoRad_PACE_OLCI.Stations.(station).Sorad; % [13x560]
    Rrs_PACE = Rrs_SoRad_PACE_OLCI.Stations.(station).PACE; % [1x14]
    Rrs_OLCI = Rrs_SoRad_PACE_OLCI.Stations.(station).OLCI; % [9x16]
    
    % Extract wavelengths for plotting
    wl_PACE = Rrs_SoRad_PACE_OLCI.wl_PACE; % [172x1]
    wl_Sorad = Rrs_SoRad_PACE_OLCI.wl_Sorad; % [560x1]
    wl_OLCI = Rrs_SoRad_PACE_OLCI.wl_OLCI; % [16x1]
    
    % Interpolate each Rrs to the common wavelength grid
    wl_common = 400:1:719; % Common wavelength grid (400 to 719 nm)
    Rrs_Sorad_interp = interp1(wl_Sorad, mean(Rrs_Sorad, 1), wl_common, 'linear', 'extrap'); % Mean across rows
    Rrs_PACE_interp = interp1(wl_PACE, Rrs_PACE, wl_common, 'linear', 'extrap');
    Rrs_OLCI_interp = interp1(wl_OLCI, mean(Rrs_OLCI, 1), wl_common, 'linear', 'extrap');
    
    % Compute standard deviation for error bars (if multiple rows for a station)
    if size(Rrs_Sorad, 1) > 1
        Rrs_Sorad_std = std(Rrs_Sorad, 0, 1); % Standard deviation across rows (along columns)
    else
        Rrs_Sorad_std = zeros(1, size(Rrs_Sorad, 2)); % No error if only 1 row
    end
    
    if size(Rrs_OLCI, 1) > 1
        Rrs_OLCI_std = std(Rrs_OLCI, 0, 1); % Standard deviation across rows (along columns)
    else
        Rrs_OLCI_std = zeros(1, size(Rrs_OLCI, 2)); % No error if only 1 row
    end
    
    % Calculate RMSE, MAPD, and Spectral Angle for the comparison
    RMSE_SoRad_PACE = sqrt(mean((Rrs_Sorad_interp - Rrs_PACE_interp).^2)); % RMSE between SoRad and PACE
    RMSE_PACE_OLCI = sqrt(mean((Rrs_PACE_interp - Rrs_OLCI_interp).^2)); % RMSE between PACE and OLCI
    MAPD_SoRad_PACE = mean(abs(Rrs_Sorad_interp - Rrs_PACE_interp) ./ Rrs_Sorad_interp) * 100; % MAPD between SoRad and PACE
    MAPD_PACE_OLCI = mean(abs(Rrs_PACE_interp - Rrs_OLCI_interp) ./ Rrs_PACE_interp) * 100; % MAPD between PACE and OLCI
    
    % Calculate Spectral Angle between SoRad and PACE
    spectral_angle_SoRad_PACE = acos(sum(Rrs_Sorad_interp .* Rrs_PACE_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_PACE_interp))) * (180 / pi); 
    
    % Spectral Angle between PACE and OLCI
    spectral_angle_PACE_OLCI = acos(sum(Rrs_PACE_interp .* Rrs_OLCI_interp) / (norm(Rrs_PACE_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);

    % Spectral Angle between SoRad and OLCI
    spectral_angle_SoRad_OLCI = acos(sum(Rrs_Sorad_interp .* Rrs_OLCI_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);
    
    % Store Spectral Angle results for the bar plot
    spectral_angles(i, :) = [spectral_angle_SoRad_PACE, spectral_angle_PACE_OLCI, spectral_angle_SoRad_OLCI];
    
    % Store RMSE and MAPD values
    RMSE_values(i, :) = [RMSE_SoRad_PACE, RMSE_PACE_OLCI];
    MAPD_values(i, :) = [MAPD_SoRad_PACE, MAPD_PACE_OLCI];
end

% Save the RMSE, MAPD, and Spectral Angles as CSV files
writetable(array2table(spectral_angles, 'VariableNames', {'SoRad_vs_PACE', 'PACE_vs_OLCI', 'SoRad_vs_OLCI'}), 'spectral_angles.csv');
writetable(array2table(RMSE_values, 'VariableNames', {'SoRad_vs_PACE', 'PACE_vs_OLCI'}), 'RMSE_values.csv');
writetable(array2table(MAPD_values, 'VariableNames', {'SoRad_vs_PACE', 'PACE_vs_OLCI'}), 'MAPD_values.csv');

%%
% Load the .mat file
load('updated_Rrs_SoRad_PACE_OLCI.mat');

% Get the station names
stations = fieldnames(Rrs_SoRad_PACE_OLCI.Stations);

% Prepare for plots
num_stations = numel(stations);

% Initialize arrays for RMSE, MAPD, and Spectral Angles
spectral_angles = zeros(num_stations, 3);
RMSE_values = zeros(num_stations, length(400:719)); % RMSE for each wavelength
MAPD_values = zeros(num_stations, length(400:719)); % MAPD for each wavelength

% Loop through each station
for i = 1:num_stations
    station = stations{i};
    
    % Extract Rrs data and wavelengths for the station
    Rrs_Sorad = Rrs_SoRad_PACE_OLCI.Stations.(station).Sorad; % [13x560]
    Rrs_PACE = Rrs_SoRad_PACE_OLCI.Stations.(station).PACE; % [1x14]
    Rrs_OLCI = Rrs_SoRad_PACE_OLCI.Stations.(station).OLCI; % [9x16]
    
    % Extract wavelengths for plotting
    wl_PACE = Rrs_SoRad_PACE_OLCI.wl_PACE; % [172x1]
    wl_Sorad = Rrs_SoRad_PACE_OLCI.wl_Sorad; % [560x1]
    wl_OLCI = Rrs_SoRad_PACE_OLCI.wl_OLCI; % [16x1]
    
    % Interpolate each Rrs to the common wavelength grid (400-719 nm)
    wl_common = 400:1:719; % Common wavelength grid (400 to 719 nm)
    Rrs_Sorad_interp = interp1(wl_Sorad, mean(Rrs_Sorad, 1), wl_common, 'linear', 'extrap'); % Mean across rows
    Rrs_PACE_interp = interp1(wl_PACE, Rrs_PACE, wl_common, 'linear', 'extrap');
    Rrs_OLCI_interp = interp1(wl_OLCI, mean(Rrs_OLCI, 1), wl_common, 'linear', 'extrap');
    
    % Compute RMSE and MAPD for each wavelength
    for j = 1:length(wl_common)
        % Get the values for the current wavelength
        SoRad_val = Rrs_Sorad_interp(j);
        PACE_val = Rrs_PACE_interp(j);
        OLCI_val = Rrs_OLCI_interp(j);
        
        % Compute RMSE for the current wavelength (SoRad vs PACE, PACE vs OLCI)
        RMSE_SoRad_PACE(j) = sqrt((SoRad_val - PACE_val)^2); % RMSE between SoRad and PACE
        RMSE_PACE_OLCI(j) = sqrt((PACE_val - OLCI_val)^2); % RMSE between PACE and OLCI
        
        % Compute MAPD for the current wavelength (SoRad vs PACE, PACE vs OLCI)
        MAPD_SoRad_PACE(j) = abs(SoRad_val - PACE_val) / SoRad_val * 100; % MAPD between SoRad and PACE
        MAPD_PACE_OLCI(j) = abs(PACE_val - OLCI_val) / PACE_val * 100; % MAPD between PACE and OLCI
    end
    
    % Store RMSE, MAPD, and Spectral Angle values for the station
    RMSE_values(i, :) = RMSE_SoRad_PACE + RMSE_PACE_OLCI; % Combine RMSE for both pairs
    MAPD_values(i, :) = MAPD_SoRad_PACE + MAPD_PACE_OLCI; % Combine MAPD for both pairs
    
    % Calculate Spectral Angle between SoRad and PACE
    spectral_angle_SoRad_PACE = acos(sum(Rrs_Sorad_interp .* Rrs_PACE_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_PACE_interp))) * (180 / pi); 
    
    % Spectral Angle between PACE and OLCI
    spectral_angle_PACE_OLCI = acos(sum(Rrs_PACE_interp .* Rrs_OLCI_interp) / (norm(Rrs_PACE_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);

    % Spectral Angle between SoRad and OLCI
    spectral_angle_SoRad_OLCI = acos(sum(Rrs_Sorad_interp .* Rrs_OLCI_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);
    
    % Store Spectral Angle results for the bar plot
    spectral_angles(i, :) = [spectral_angle_SoRad_PACE, spectral_angle_PACE_OLCI, spectral_angle_SoRad_OLCI];
end

% Save the RMSE, MAPD, and Spectral Angles as CSV files
%%
% Load the .mat file
load('updated_Rrs_SoRad_PACE_OLCI.mat');

% Get the station names
stations = fieldnames(Rrs_SoRad_PACE_OLCI.Stations);

% Prepare for plots
num_stations = numel(stations);

% Initialize arrays for spectral angles, RMSE, and MAPD
spectral_angles = zeros(num_stations, 3);
RMSE_values = zeros(num_stations, length(400:719)); % RMSE for each wavelength
MAPD_values = zeros(num_stations, length(400:719)); % MAPD for each wavelength

% Define the wavelength range
wl_common = 400:1:719; % Common wavelength grid (400 to 719 nm)

% Loop through each station
for i = 1:num_stations
    station = stations{i};
    
    % Extract Rrs data and wavelengths for the station
    Rrs_Sorad = Rrs_SoRad_PACE_OLCI.Stations.(station).Sorad; % [13x560]
    Rrs_PACE = Rrs_SoRad_PACE_OLCI.Stations.(station).PACE; % [1x14]
    Rrs_OLCI = Rrs_SoRad_PACE_OLCI.Stations.(station).OLCI; % [9x16]
    
    % Extract wavelengths for plotting
    wl_PACE = Rrs_SoRad_PACE_OLCI.wl_PACE; % [172x1]
    wl_Sorad = Rrs_SoRad_PACE_OLCI.wl_Sorad; % [560x1]
    wl_OLCI = Rrs_SoRad_PACE_OLCI.wl_OLCI; % [16x1]
    
    % Interpolate each Rrs to the common wavelength grid (400-719 nm)
    Rrs_Sorad_interp = interp1(wl_Sorad, mean(Rrs_Sorad, 1), wl_common, 'linear', 'extrap'); % Mean across rows
    Rrs_PACE_interp = interp1(wl_PACE, Rrs_PACE, wl_common, 'linear', 'extrap');
    Rrs_OLCI_interp = interp1(wl_OLCI, mean(Rrs_OLCI, 1), wl_common, 'linear', 'extrap');
    
    % Compute RMSE and MAPD for each wavelength
    for j = 1:length(wl_common)
        % Get the values for the current wavelength
        SoRad_val = Rrs_Sorad_interp(j);
        PACE_val = Rrs_PACE_interp(j);
        OLCI_val = Rrs_OLCI_interp(j);
        
        % Compute RMSE for the current wavelength (SoRad vs PACE, PACE vs OLCI)
        RMSE_SoRad_PACE(j) = sqrt((SoRad_val - PACE_val)^2); % RMSE between SoRad and PACE
        RMSE_PACE_OLCI(j) = sqrt((PACE_val - OLCI_val)^2); % RMSE between PACE and OLCI
        
        % Compute MAPD for the current wavelength (SoRad vs PACE, PACE vs OLCI)
        MAPD_SoRad_PACE(j) = abs(SoRad_val - PACE_val) / SoRad_val * 100; % MAPD between SoRad and PACE
        MAPD_PACE_OLCI(j) = abs(PACE_val - OLCI_val) / PACE_val * 100; % MAPD between PACE and OLCI
    end
    
    % Store RMSE, MAPD, and Spectral Angle values for the station
    RMSE_values(i, :) = RMSE_SoRad_PACE + RMSE_PACE_OLCI; % Combine RMSE for both pairs
    MAPD_values(i, :) = MAPD_SoRad_PACE + MAPD_PACE_OLCI; % Combine MAPD for both pairs
    
    % Calculate Spectral Angle between SoRad and PACE
    spectral_angle_SoRad_PACE = acos(sum(Rrs_Sorad_interp .* Rrs_PACE_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_PACE_interp))) * (180 / pi); 
    
    % Spectral Angle between PACE and OLCI
    spectral_angle_PACE_OLCI = acos(sum(Rrs_PACE_interp .* Rrs_OLCI_interp) / (norm(Rrs_PACE_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);

    % Spectral Angle between SoRad and OLCI
    spectral_angle_SoRad_OLCI = acos(sum(Rrs_Sorad_interp .* Rrs_OLCI_interp) / (norm(Rrs_Sorad_interp) * norm(Rrs_OLCI_interp))) * (180 / pi);
    
    % Store Spectral Angle results for the bar plot
    spectral_angles(i, :) = [spectral_angle_SoRad_PACE, spectral_angle_PACE_OLCI, spectral_angle_SoRad_OLCI];
end

% Convert wavelength values (numeric) to cell array of strings for the table column names
wavelength_strings = cellfun(@num2str, num2cell(wl_common), 'UniformOutput', false);

% Save the RMSE, MAPD, and Spectral Angles as CSV files
writetable(array2table(spectral_angles, 'VariableNames', {'SoRad_vs_PACE', 'PACE_vs_OLCI', 'SoRad_vs_OLCI'}), 'spectral_angles.csv');
writetable(array2table(RMSE_values, 'VariableNames', wavelength_strings), 'RMSE_values.csv');
writetable(array2table(MAPD_values, 'VariableNames', wavelength_strings), 'MAPD_values.csv');
