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
exportgraphics(gcf,filename,'Resolution',400)