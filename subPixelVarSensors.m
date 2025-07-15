% Clear workspace and command window
clear;
clc;

% Load the .mat file
matFilePath = '/Users/masud/OneDriveUGA/CruiseRVSavannah/Subpixel/updated_Rrs_SoRad_PACE_OLCI.mat';
if ~isfile(matFilePath)
    error('File %s not found', matFilePath);
end
data = load(matFilePath);
if ~isfield(data, 'Rrs_SoRad_PACE_OLCI')
    error('Rrs_SoRad_PACE_OLCI structure not found in %s', matFilePath);
end
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
    
    % Plot SoRad with error bars (if multiple rows)
    if size(stData.Sorad, 1) > 1
        errorbar(RrsStruct.wl_Sorad, soradMean, soradStd, 'b-', 'DisplayName', 'SoRad');
    else
        plot(RrsStruct.wl_Sorad, soradMean, 'b-', 'DisplayName', 'SoRad');
    end
    
    % Interpolate PACE to SoRad wavelengths (400-800 nm)
    idx_wl = RrsStruct.wl_Sorad >= 400 & RrsStruct.wl_Sorad <= 800;
    wl_common = RrsStruct.wl_Sorad(idx_wl);
    paceMean_interp = interp1(RrsStruct.wl_PACE, paceMean, wl_common, 'linear', 'extrap');
    paceStd_interp = interp1(RrsStruct.wl_PACE, paceStd, wl_common, 'linear', 'extrap');
    if size(stData.PACE, 1) > 1
        errorbar(wl_common, paceMean_interp, paceStd_interp, 'r-', 'DisplayName', 'PACE');
    else
        plot(wl_common, paceMean_interp, 'r-', 'DisplayName', 'PACE');
    end
    
    % Interpolate OLCI to SoRad wavelengths
    olciMean_interp = interp1(RrsStruct.wl_OLCI, olciMean, wl_common, 'linear', 'extrap');
    olciStd_interp = interp1(RrsStruct.wl_OLCI, olciStd, wl_common, 'linear', 'extrap');
    if size(stData.OLCI, 1) > 1
        errorbar(wl_common, olciMean_interp, olciStd_interp, 'g-', 'DisplayName', 'OLCI');
    else
        plot(wl_common, olciMean_interp, 'g-', 'DisplayName', 'OLCI');
    end
    
    % Plot settings
    xlabel('Wavelength (nm)');
    ylabel('Rrs (sr^{-1})');
    title(sprintf('Station %s', stName));
    legend('show');
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

% Bar plot for Spectral Angles
figure('Name', 'Spectral Angles by Station', 'Position', [100, 100, 800, 400]);
barData = [[metrics.SoRad_vs_PACE_SA]; [metrics.SoRad_vs_OLCI_SA]; [metrics.PACE_vs_OLCI_SA]]';
bar(stations, barData);
xlabel('Station');
ylabel('Spectral Angle (degrees)');
title('Spectral Angles Between Sensors');
legend('SoRad vs PACE', 'SoRad vs OLCI', 'PACE vs OLCI', 'Location', 'best');
grid on;