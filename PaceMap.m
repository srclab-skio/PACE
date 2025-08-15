clear all;
close all;
clc;
% Specify the input NetCDF file
% input_file = '/Volumes/SRC_HDD_Mas/Data/PACE_Sapelo_Matchups/PACE_Sapelo_Matchups/PACE_OCI.20250726T181220.L2.OC_BGC.V3_0.NRT_subset.nc';
input_file = '/Volumes/SRC_HDD_Mas/Data/PACE_Sapelo_Matchups/PACE_Sapelo_Matchups/PACE_OCI.20250726T181220.L2.OC_AOP.V3_0.NRT_subset.nc';
% ncdisp("PACE_OCI.20250726T181220.L2.OC_AOP.V3_0.NRT_subset.nc");

% Read in latitude, longitude, and the variable you want to plot (e.g., 'ocn_data')
lat = ncread(input_file, 'lat');  % Replace 'lat' with the correct latitude variable name
lon = ncread(input_file, 'lon');  % Replace 'lon' with the correct longitude variable name
% data_var = ncread(input_file, 'chlor_a');  % Replace 'your_variable' with the variable to plot
data_var = ncread(input_file, 'Rrs_670');
% Check the sizes of lat, lon, and data_var
% size(lat)
% size(lon)
% size(data_var)

% Flatten lat and lon into 1D vectors (if necessary)
% Use reshape to convert the 2D lat and lon grids into 1D vectors
% lat_vec = reshape(lat, [], 1);
% lon_vec = reshape(lon, [], 1);
% data_var_vec = reshape(data_var, [], 1);
% Replace NaN values with a specific value (0, which corresponds to black in the colormap)
% data_var(isnan(data_var)) = 0;
% Make sure data_var is also reshaped to match the lat-lon grid

%%
figure;
 m_proj('lambert','lat',[28 34],'lon',[-78. -81.5]);
 m_pcolor(lon, lat, data_var);
 m_grid('tickdir','out', ...
        'linestyle','none','ticklen',.02,'FontName',"Leelawadee", 'bold');
%  m_grid('tickdir','out','ytick',[ 31.5 32 33 34],'xtick',[-81 -80 -79 -78 -77.5], ...
%         'linestyle','none','ticklen',.02,'FontName',"Leelawadee", 'bold');
 m_gshhs_h('patch',[.5 .5 .5]);

 colormap(jet)
 set(gca,'ColorScale','log')
%  clim([-1 5])

 m_text(-80.9,33.87,"t_0: " ,'FontSize',7,'FontName',"Leelawadee")
 m_text(-80.9,33.67,"t_1: " ,'FontSize',7,'FontName',"Leelawadee")
 title("PACE OCI",'FontSize',12,'FontName',"Leelawadee")
 set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold');


c = colorbar('southoutside','FontSize',12);
% c.Position(1) = c.Position(1)-0.435; % Why this position here, and th enumber -0.435
% c.Position(2) = c.Position(2)-0.10;
% c.Position(3) = c.Position(3)*2.5;
% c.Position(4) = c.Position(4)*1.5;
set(c,'FontSize', 12,'FontName','Leelawadee')
set(c,'ytick',([0.2 0.5 1 2 5]),'yticklabel',["<0.2","0.5", "1", "2", "5"],'tickdir','out')
set(get(c,'xlabel'),'string','Log-Chlorophyll-a [mg/m^3]', 'FontSize', 12,'FontName','Leelawadee')
% figNames  = ("SAB_Chlor_a_" + band_head(1) + S.SH.t0(1:10) + ".png");
% Set font weight and size for axis labels and axis line thickness
%%%%%%%%%%%
set(get(gca, 'XLabel'), 'FontSize', 12, 'FontWeight', 'bold');
set(get(gca, 'YLabel'), 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2);
set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold');

filename = fullfile('/Users/masud/Documents/MATLAB/SeaHawk_comparison/surfPlotsAll',"Chla" + satLocation  + S.SH.t0(1:10) + ".png");
% export figure
% exportgraphics(gcf,("AG_Chlor_a_" + band_head + S.SH.t0(1:10) + ".png"),'Resolution',600)
% exportgraphics(gcf,'myfigure.eps','ContentType',"vector")
% The next 2 lines should be uncommented for Chla and saving figs
% exportgraphics(gcf,filenameEps,'Resolution',400)
% exportgraphics(gcf,filename,'Resolution',400)
<<<<<<< Updated upstream
exportgraphics(gcf,filename,'Resolution',400)
=======
% exportgraphics(gcf,filename,'Resolution',400)


%% TSS MAP
clear all;
% Process PACE NetCDF data to calculate and plot Total Suspended Solids (TSS)
% Process PACE NetCDF data to calculate and plot Total Suspended Solids (TSS)
% Step 1: Define file path
nc_file = '/Users/masud/OneDriveUGA/QWIP/subset_PACE.nc';  % Updated to full path

ncinfo(nc_file)
% Specify the NetCDF file path
% Get info about the NetCDF file
nc_info = ncinfo(nc_file);

% Extract all variable names
all_vars = {nc_info.Variables.Name};

% Find variables starting with 'Rrs_'
rrs_vars = all_vars(startsWith(all_vars, 'Rrs_'));

% Display the result
disp(rrs_vars);
bbox = [29, 32, -81, -79];


%% % MATLAB script to calculate TSS from Rrs_665 in the NC file and visualize it using m_map

% Note: This script assumes that the m_map toolbox is installed and available in your MATLAB path.
% You can download m_map from: https://www.eoas.ubc.ca/~rich/map.html

% Define the input file path
input_file = '/Volumes/SRC_HDD_Mas/Data/PACE_Sapelo_Matchups/subset_PACE.nc';
% Get info about the NetCDF file
nc_info = ncinfo(input_file);
% Extract all variable names
all_vars = {nc_info.Variables.Name};
% Find variables starting with 'Rrs_'
rrs_vars = all_vars(startsWith(all_vars, 'Rrs_'));

% Display the result
disp(rrs_vars);

% Read the raw Rrs_665 data (int16)
rrs_raw = ncread(input_file, 'Rrs_555');

% Define scaling parameters from attributes
scale_factor = 2e-6;
add_offset = 0.05;
fill_value = -32767;

% Apply scaling and offset to get actual Rrs values
% rrs = double(rrs_raw) * scale_factor + add_offset;
rrs = double(rrs_raw) ;
% Set fill values to NaN
rrs(rrs_raw == fill_value) = NaN;

% Calculate TSS using the provided algorithm
% tss = 1.74 + (355.85 * pi * rrs) ./ (1 - (pi * rrs) / 1728);
tss = 2 + (500 * pi * rrs) ./ (1 - (pi * rrs) / 1728);


% Read latitude and longitude grids
lat = ncread(input_file, 'lat');
lon = ncread(input_file, 'lon');

% Compute map limits
% lonmin = min(lon(:));
% lonmax = max(lon(:));
% latmin = min(lat(:));
% latmax = max(lat(:));
lonmin = -82;
lonmax = -78;
latmin = 29;
latmax = 34;

% Visualize the TSS map using m_map with Lambert projection
figure;
m_proj('lambert', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);

% Plot the TSS data
m_pcolor(lon, lat, tss);
shading flat;

% Add high-resolution coastline, gray land, and internal waters (lakes as holes)
hold on;
m_gshhs_h('patch', [0.7 0.7 0.7], 'edgecolor', 'k');  % Gray land fill with black coastline

% Add grid and labels (optional for better visualization)
m_grid('box', 'fancy', 'tickdir', 'in');

% Colorbar and title
colorbar;
title('Total Suspended Solids (TSS) in mg L^{-1}');
caxis([min(tss(:), [], 'omitnan') max(tss(:), [], 'omitnan')]);  % Fit colorbar to TSS value range
colormap jet;  % Optional: Choose a colormap, e.g., jet, parula, etc.

%%
% MATLAB script to calculate TSS from Rrs_488, Rrs_555, and Rrs_658 in the NC file and visualize it using m_map

% Note: This script assumes that the m_map toolbox is installed and available in your MATLAB path.
% You can download m_map from: https://www.eoas.ubc.ca/~rich/map.html
% Using Rrs_658 (658.34 nm) as a substitute for Rrs_645, as Rrs_645 is not available in the provided NetCDF file.

% Define the input file path
input_file = '/Volumes/SRC_HDD_Mas/Data/PACE_Sapelo_Matchups/subset_PACE_Rrs.nc';

% Read the raw Rrs data (int16) for required wavelengths
rrs_488_raw = ncread(input_file, 'Rrs_487');  % Closest to 488 nm
rrs_555_raw = ncread(input_file, 'Rrs_555');
rrs_645_raw = ncread(input_file, 'Rrs_658');  % Closest to 645 nm

% Define scaling parameters from attributes
scale_factor = 2e-6;
add_offset = 0.05;
fill_value = -32767;

% Apply scaling and offset to get actual Rrs values
% rrs_488 = double(rrs_488_raw) * scale_factor + add_offset;
% rrs_555 = double(rrs_555_raw) * scale_factor + add_offset;
% rrs_645 = double(rrs_645_raw) * scale_factor + add_offset;

rrs_488 = double(rrs_488_raw) ;
rrs_555 = double(rrs_555_raw) ;
rrs_645 = double(rrs_645_raw) ;
% Set fill values to NaN
rrs_488(rrs_488_raw == fill_value) = NaN;
rrs_555(rrs_555_raw == fill_value) = NaN;
rrs_645(rrs_645_raw == fill_value) = NaN;

% Calculate log10(TSS) using the provided equation
tss = 0.6311 + 22.2158 * (rrs_555 + rrs_645) - 0.5239 * (rrs_488 ./ rrs_555);

% Convert to TSS
tss = 10 .^ log10TSS;

% Read latitude and longitude grids
lat = ncread(input_file, 'lat');
lon = ncread(input_file, 'lon');

% Compute map limits
lonmin = min(lon(:));
lonmax = max(lon(:));
latmin = min(lat(:));
latmax = max(lat(:));

% Visualize the TSS map using m_map with Lambert projection
figure;
m_proj('lambert', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);

% Plot the TSS data
m_pcolor(lon, lat, tss);
shading flat;

% Add high-resolution coastline, gray land, and internal waters (lakes as holes)
hold on;
m_gshhs_h('patch', [0.7 0.7 0.7], 'edgecolor', 'k');  % Gray land fill with black coastline

% Add grid and labels
m_grid('box', 'fancy', 'tickdir', 'in');

% Colorbar and title
colorbar;
title('Total Suspended Solids (TSS) in mg L^{-1}');
caxis([min(tss(:), [], 'omitnan') max(tss(:), [], 'omitnan')]);  % Fit colorbar to TSS value range
colormap jet;  % Optional: Choose a colormap, e.g., jet, parula, etc.


>>>>>>> Stashed changes
