clear all; close all; clc;
% Load the .mat file (adjust the file path)
load('/Users/masud/OneDriveUGA/CruiseRVSavannah/Comparison_plots/SVC_PACE_Rrs_7July2025.mat');  % Replace with the actual path to your .mat file

% Read the wavelengths and corresponding Rrs values for each spectrum
wavelength1 = Rrs_final_SVC_PACE.wl;  % Wavelengths for Rrs_3C
Rrs1 = Rrs_final_SVC_PACE.Rrs_3C;  % Rrs values for Rrs_3C

wavelength2 = Rrs_final_SVC_PACE.wl_PACE;  % Wavelengths for Rrs_PACE
Rrs2 = Rrs_final_SVC_PACE.Rrs_PACE;  % Rrs values for Rrs_PACE

wavelength3 = Rrs_final_SVC_PACE.wl_Sorad;  % Wavelengths for Rrs_Sorad
Rrs3 = Rrs_final_SVC_PACE.Rrs_Sorad;  % Rrs values for Rrs_Sorad

wavelength4 = Rrs_final_SVC_PACE.wl_HPro;  % Wavelengths for Rrs_HPro
Rrs4 = Rrs_final_SVC_PACE.Rrs_HPro;  % Rrs values for Rrs_HPro
wavelength5 = Rrs_final_SVC_PACE.wl_OLCI;  % Wavelengths for Rrs_HPro
Rrs5 = Rrs_final_SVC_PACE.Rrs_OLCI;  % Rrs values for Rrs_HPro

% Define the common wavelength range (350 nm to 800 nm)
wavelength_range = 350:800;  % Wavelength range for the x-axis


%% For each station
% Create a figure for 16 subplots
figure;
% set(gcf, 'Position', [100, 100, 800, 800]);  % Adjust the figure size
% Ensure saved PDF matches screen size
set(gcf, 'PaperPositionMode', 'auto'); 
sgtitle('Comparison for all sensors within 2hr of PACE overpass', 'FontSize', 16, 'FontWeight', 'bold');  % Set figure title

% Loop through each of the 16 sites (rows)
for i = [1:16]
    % Create a subplot for the current sampling site
    %subplot(4, 4, i);  % 4x4 grid of subplots (16 total)
    subplot(4, 4, i);  % 4x4 grid of subplots (16 total)

    % Plot the spectra for the current site (i-th row) with different colors
    plot(wavelength1, Rrs1(i,:), 'r-', 'LineWidth', 1.5); % Red for Rrs_3C
    hold on;
    plot(wavelength2, Rrs2(i,:), 'g-', 'LineWidth', 1.5); % Green for Rrs_PACE
    plot(wavelength3, Rrs3(i,:), 'b-', 'LineWidth', 1.5); % Blue for Rrs_Sorad
    plot(wavelength4, Rrs4(i,:), 'k-', 'LineWidth', 1.5); % Black for Rrs_HPro
    % plot(wavelength5, Rrs5(i,:), 'c-', 'LineWidth', 1.5); % Black for Rrs_OLCI
    plot(wavelength5, Rrs5(i,:), 'c-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k'); 

    % Customize each subplot
    xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Rrs (sr^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([350, 800]);  % Set x-axis limits from 350 nm to 800 nm
    ylim([-0.001, 0.012]);  % Set x-axis limits from 350 nm to 800 nm
    title(['Station ' num2str(i)], 'FontSize', 14, 'FontWeight', 'bold');  % Title for each subplot
    grid on;  % Enable grid
    set(gca, 'GridLineStyle', '-');  % Set grid lines
    set(gca, 'LineWidth', 1.5);  % Increase axis line width
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');  % Make tick labels bold and larger
    % Remove the grid (graticule) from the figure
    grid off;  % Turn off the grid

    % Add legend to only the first subplot (i == 1)
    if i == 16
        legend('SVC', 'PACE', 'SoRad', 'HPro','OLCI-Acolite', 'Location', 'northeast', 'FontSize', 10, 'FontWeight', 'bold');
    end

%     % Add legend to each subplot
%     legend('SVC', 'PACE', 'Sorad', 'HPro', 'Location', 'Best', 'FontSize', 10, 'FontWeight', 'bold');
end
    % Add legend to each subplot
    %legend('SVC', 'PACE', 'Sorad', 'HPro', 'Location', 'northeast', 'FontSize', 10, 'FontWeight', 'bold');
% Adjust figure appearance
set(gca, 'TickDir', 'out');  % Place ticks outside the plot
set(gca, 'TickLength', [0.02, 0.02]);  % Adjust tick length for better visibility

% Save the figure with high resolution (400 dpi) and vector image (PDF)
% saveas(gcf, 'overlayed_spectra_400dpi_OLCI.png');  % Save as PNG with high resolution
% saveas(gcf, 'overlayed_spectra_OLCI.eps');  % Save as PDF (vector image)
 % exportgraphics(gcf, 'all_Sts_all_Sensors.pdf', 'ContentType','vector')

%% FOr only good data with PACE window
clear all; close all; clc;

% Load the .mat file (adjust the file path)
% load('/Users/masud/OneDriveUGA/CruiseRVSavannah/Comparison_plots/SVC_PACE_Rrs_7July2025.mat');  % Replace with the actual path to your .mat file
% Load the .mat file (adjust the file path)
load('/Users/masud/OneDriveUGA/CruiseRVSavannah/Comparison_plots/SVC_PACE_Rrs_4Aug2025.mat');  % Replace with the actual path to your .mat file

% Read the wavelengths and corresponding Rrs values for each spectrum
wavelength1 = Rrs_final_SVC_PACE.wl;  % Wavelengths for Rrs_3C
Rrs1 = Rrs_final_SVC_PACE.Rrs_3C;  % Rrs values for Rrs_3C

wavelength2 = Rrs_final_SVC_PACE.wl_PACE;  % Wavelengths for Rrs_PACE
Rrs2 = Rrs_final_SVC_PACE.Rrs_PACE;  % Rrs values for Rrs_PACE

wavelength3 = Rrs_final_SVC_PACE.wl_Sorad;  % Wavelengths for Rrs_Sorad
Rrs3 = Rrs_final_SVC_PACE.Rrs_Sorad;  % Rrs values for Rrs_Sorad

wavelength4 = Rrs_final_SVC_PACE.wl_HPro;  % Wavelengths for Rrs_HPro
Rrs4 = Rrs_final_SVC_PACE.Rrs_HPro;  % Rrs values for Rrs_HPro
wavelength5 = Rrs_final_SVC_PACE.wl_OLCI;  % Wavelengths for Rrs_HPro
Rrs5 = Rrs_final_SVC_PACE.Rrs_OLCI_Aco;  % Rrs values for Rrs_HPro

% Define the common wavelength range (350 nm to 800 nm)
wavelength_range = 350:800;  % Wavelength range for the x-axis

%% For all station on 30th
% Create a figure for 6 subplots (stations 10 to 16 excluding station 13)
figure;
% set(gcf, 'Position', [100, 100, 800, 800]);  % Adjust the figure size
% Ensure saved PDF matches screen size
set(gcf, 'PaperPositionMode', 'auto'); 
sgtitle('Comparison for stations 10 to 16', 'FontSize', 16, 'FontWeight', 'bold');  % Set figure title

% Loop through stations 10 to 16, excluding station 13
station_indices = [10:13, 15:16];  % Exclude station 13# Currently St 14, cuz 9.5 is st 10
for i = 1:length(station_indices)
    % Create a subplot for the current sampling site
    subplot(2, 3, i);  % 2x3 grid of subplots (6 total)

    % Plot the spectra for the current site (i-th row) with different colors
    plot(wavelength1, Rrs1(station_indices(i), :), 'Color', [0 0 0.5], 'LineWidth', 1.5); % Red for Rrs_3C
    hold on;
    plot(wavelength2, Rrs2(station_indices(i), :), 'g-', 'LineWidth', 1.5); % Green for Rrs_PACE
    plot(wavelength3, Rrs3(station_indices(i), :), 'b-', 'LineWidth', 1.5); % Blue for Rrs_Sorad
    plot(wavelength4, Rrs4(station_indices(i), :), 'k-', 'LineWidth', 1.5); % Black for Rrs_HPro
    plot(wavelength5, Rrs5(station_indices(i), :), 'c-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k'); % Cyan for Rrs_OLCI

    % Customize each subplot
    xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Rrs (sr^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([350, 800]);  % Set x-axis limits from 350 nm to 800 nm
    ylim([-0.001, 0.009]);  % Set y-axis limits
    title(['St. ' num2str(station_indices(i))], 'FontSize', 14, 'FontWeight', 'bold');  % Title for each subplot
    grid off;  % Enable grid
    set(gca, 'GridLineStyle', '-');  % Set grid lines
    set(gca, 'LineWidth', 1.5);  % Increase axis line width
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');  % Make tick labels bold and larger

    % Add legend to only the last subplot
    if i == length(station_indices)
        legend('SVC', 'PACE', 'SoRad', 'HPro', 'OLCI-Acolite', 'Location', 'northeast', 'FontSize', 10, 'FontWeight', 'bold');
    end
end

% Adjust figure appearance
set(gca, 'TickDir', 'out');  % Place ticks outside the plot
set(gca, 'TickLength', [0.02, 0.02]);  % Adjust tick length for better visibility

% Save the figure with high resolution (400 dpi) and vector image (PDF)
% saveas(gcf, 'overlayed_spectra_400dpi_OLCI.png');  % Save as PNG with high resolution
% saveas(gcf, 'overlayed_spectra_OLCI.eps');  % Save as PDF (vector image)
% exportgraphics(gcf, 'stations_10_16_excluding_13.pdf', 'ContentType', 'vector');
filename = fullfile('/Users/masud/OneDriveUGA/CruiseRVSavannah/Comparison_plots/', "All_Good_Pace_OLCI_Acolite" + ""  + ".png");
exportgraphics(gcf,filename,'Resolution',600)


%%% For only flagged and valid Data, St 12, 14, 15
%
%clear all; close all; clc;
%
%% Load the .mat file (adjust the file path)
%load('/Users/masud/OneDriveUGA/CruiseRVSavannah/Comparison_plots/SVC_PACE_Rrs_7July2025.mat');  % Replace with the actual path to your .mat file
%
%% Read the wavelengths and corresponding Rrs values for each spectrum
%wavelength1 = Rrs_final_SVC_PACE.wl;  % Wavelengths for Rrs_3C
%Rrs1 = Rrs_final_SVC_PACE.Rrs_3C;  % Rrs values for Rrs_3C
%
%wavelength2 = Rrs_final_SVC_PACE.wl_PACE;  % Wavelengths for Rrs_PACE
%Rrs2 = Rrs_final_SVC_PACE.Rrs_PACE;  % Rrs values for Rrs_PACE
%
%wavelength3 = Rrs_final_SVC_PACE.wl_Sorad;  % Wavelengths for Rrs_Sorad
%Rrs3 = Rrs_final_SVC_PACE.Rrs_Sorad;  % Rrs values for Rrs_Sorad
%
%wavelength4 = Rrs_final_SVC_PACE.wl_HPro;  % Wavelengths for Rrs_HPro
%Rrs4 = Rrs_final_SVC_PACE.Rrs_HPro;  % Rrs values for Rrs_HPro
%wavelength5 = Rrs_final_SVC_PACE.wl_OLCI;  % Wavelengths for Rrs_HPro
%Rrs5 = Rrs_final_SVC_PACE.Rrs_OLCI;  % Rrs values for Rrs_HPro
%
%% Define the common wavelength range (350 nm to 800 nm)
%wavelength_range = 350:800;  % Wavelength range for the x-axis
%
%%% Create a figure for 3 subplots (stations 12, 14, and 15)
%figure;
%
%% set(gcf, 'Position', [100, 100, 800, 800]);  % Adjust the figure size
%% Ensure saved PDF matches screen size
%set(gcf, 'PaperPositionMode', 'auto'); 
%sgtitle('Comparison for stations 12, 14, and 15', 'FontSize', 16, 'FontWeight', 'bold');  % Set figure title
%
%% Loop through stations 12, 14, and 15
%station_indices = [12, 14, 15];  % Stations 12, 14, and 15
%for i = 1:length(station_indices)
%    % Create a subplot for the current sampling site
%    subplot(1, 3, i);  % 1x3 grid of subplots (3 total)
%
%    % Plot the spectra for the current site (i-th row) with different colors
%    plot(wavelength1, Rrs1(station_indices(i), :), 'r-', 'LineWidth', 1.5); % Red for Rrs_3C
%    hold on;
%    plot(wavelength2, Rrs2(station_indices(i), :), 'g-', 'LineWidth', 1.5); % Green for Rrs_PACE
%    plot(wavelength3, Rrs3(station_indices(i), :), 'b-', 'LineWidth', 1.5); % Blue for Rrs_Sorad
%    plot(wavelength4, Rrs4(station_indices(i), :), 'k-', 'LineWidth', 1.5); % Black for Rrs_HPro
%    plot(wavelength5, Rrs5(station_indices(i), :), 'c-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k'); % Cyan for Rrs_OLCI
%
%    % Customize each subplot
%    xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight', 'bold');
%    ylabel('Rrs (sr^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
%    xlim([350, 800]);  % Set x-axis limits from 350 nm to 800 nm
%    ylim([-0.001, 0.009]);  % Set y-axis limits
%    title(['St. ' num2str(station_indices(i))], 'FontSize', 14, 'FontWeight', 'bold');  % Title for each subplot
%    grid off;  % Enable grid
%    set(gca, 'GridLineStyle', '-');  % Set grid lines
%    set(gca, 'LineWidth', 1.5);  % Increase axis line width
%    set(gca, 'FontSize', 12, 'FontWeight', 'bold');  % Make tick labels bold and larger
%
%    % Add legend to only the last subplot
%    if i == length(station_indices)
%        legend('SVC', 'PACE', 'SoRad', 'HPro', 'OLCI', 'Location', 'northeast', 'FontSize', 10, 'FontWeight', 'bold');
%    end
%end
%
%% Adjust figure appearance
%set(gca, 'TickDir', 'out');  % Place ticks outside the plot
%set(gca, 'TickLength', [0.02, 0.02]);  % Adjust tick length for better visibility
%% Ensure figure size fits the subplots properly (auto-adjusted for space)
%set(gcf, 'Position', [100, 100, 1000, 400]);  % Adjust figure size if needed
%% Automatically adjust the layout of subplots for better spacing
%% tight_layout();
%% exportgraphics(gcf, 'stations_12_14_15.pdf', 'ContentType', 'vector');
%filename = fullfile('/Users/masud/OneDriveUGA/CruiseRVSavannah/Comparison_plots/', "QC Controlled" + "3"  + ".png");
%exportgraphics(gcf,filename,'Resolution',400)















