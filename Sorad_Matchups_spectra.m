% %% Step 1: Initialize Workspace and Set Directory
% %% 21 October 2026
% clc; clear; close all;
% 
% % Set Working Directory
% cd('/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/CleanedRrs');
% 
% % Define Input File Paths
% data1_path = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/St_cruise_insitu_TSS.csv';
% data2_path = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad/rrs_with_rrswl_09_2_SAV.csv'; 
% output_csv = 'St_SRC_Optics_cruise_PACE_Matchups.csv';
% 
% %% Step 2: Load Data Files with Explicit Header Parsing Options
% fprintf('Loading Data 1 (Station In-Situ TSS): %s\n', data1_path);
% opts1 = detectImportOptions(data1_path);
% opts1.VariableNamingRule = 'preserve';
% data1 = readtable(data1_path, opts1);
% 
% fprintf('Loading Data 2 (SoRad Spectrometer Data): %s\n', data2_path);
% opts2 = detectImportOptions(data2_path);
% opts2.VariableNamingRule = 'preserve';
% opts2.VariableNamesLine = 1; % Force Row 1 as column headers
% opts2.DataLines = [2 Inf];   % Data starts on Row 2
% data2 = readtable(data2_path, opts2);
% 
% %% Step 3: Standardize SoRad Wavelength Headers (340..899 -> Rrs_340..Rrs_899)
% vars2 = data2.Properties.VariableNames;
% numeric_band_idx = [];
% new_vars2 = vars2;
% 
% for i = 1:length(vars2)
%     var_str = char(strtrim(vars2{i}));
%     % Clean MATLAB auto-prefixed characters ('x', quotes, etc.)
%     var_clean = regexprep(var_str, '[^0-9a-zA-Z_]', '');
%     tokens = regexp(var_clean, '(\d+)', 'tokens');
% 
%     if ~isempty(tokens)
%         wl_num = str2double(tokens{1}{1});
%         % Identify hyperspectral channels (340 nm to 899 nm)
%         if wl_num >= 340 && wl_num <= 899
%             new_vars2{i} = sprintf('Rrs_%d', round(wl_num));
%             numeric_band_idx = [numeric_band_idx, i]; %#ok<AGROW>
%         end
%     end
% end
% 
% % Re-assign normalized header names back to Data 2 table
% data2.Properties.VariableNames = new_vars2;
% rrs_headers_data2 = new_vars2(numeric_band_idx);
% 
% fprintf('Successfully identified %d Rrs channels in SoRad Data 2 (%s to %s).\n', ...
%     length(rrs_headers_data2), rrs_headers_data2{1}, rrs_headers_data2{end});
% 
% %% Step 4: Dynamically Match Key Metadata Column Names
% vars2_clean = data2.Properties.VariableNames;
% 
% % Find Date column
% date_col_idx = find(contains(vars2_clean, 'date', 'IgnoreCase', true), 1);
% % Find Time column
% time_col_idx = find(contains(vars2_clean, 'time', 'IgnoreCase', true) | ...
%                     contains(vars2_clean, 'utc', 'IgnoreCase', true), 1);
% % Find Latitude and Longitude columns
% lat_col_idx = find(contains(vars2_clean, 'lat', 'IgnoreCase', true), 1);
% lon_col_idx = find(contains(vars2_clean, 'lon', 'IgnoreCase', true), 1);
% 
% if isempty(date_col_idx) || isempty(time_col_idx) || isempty(lat_col_idx) || isempty(lon_col_idx)
%     error('Could not locate required metadata columns (Date, Time, Lat, Lon) in Data 2.');
% end
% 
% date_col_name = vars2_clean{date_col_idx};
% time_col_name = vars2_clean{time_col_idx};
% lat_col_name  = vars2_clean{lat_col_idx};
% lon_col_name  = vars2_clean{lon_col_idx};
% 
% fprintf('Data 2 metadata columns -> Date: "%s", Time: "%s", Lat: "%s", Lon: "%s"\n', ...
%     date_col_name, time_col_name, lat_col_name, lon_col_name);
% 
% %% Step 5: Convert and Align UTC Timestamps to EST (UTC - 5 hours)
% % Parse Data 1 EST Timestamps
% data1_dt_str = strcat(cellstr(string(data1.Date)), {' '}, cellstr(string(data1.TimeEST)));
% data1_datetime = datetime(data1_dt_str, 'InputFormat', 'MM/dd/yy H:mm', 'TimeZone', 'America/New_York');
% 
% % Parse Data 2 UTC Timestamps & Convert to EST
% data2_dt_str = strcat(cellstr(string(data2.(date_col_name))), {' '}, cellstr(string(data2.(time_col_name))));
% data2_datetime_utc = datetime(data2_dt_str, 'InputFormat', 'MM/dd/yy HH:mm:ss', 'TimeZone', 'UTC');
% data2_datetime_est = datetime(data2_datetime_utc, 'TimeZone', 'America/New_York');
% 
% %% Step 6: Spatio-Temporal Matchup & Ensemble Averaging Loop
% N_stations = height(data1);
% N_rrs_bands = length(rrs_headers_data2);
% matched_rrs_matrix = NaN(N_stations, N_rrs_bands);
% data2_rrs_array = table2array(data2(:, rrs_headers_data2));
% 
% % Extract coordinates as double arrays
% data2_lats = double(data2.(lat_col_name));
% data2_lons = double(data2.(lon_col_name));
% 
% % Matchup Tolerances
% geo_tolerance_deg = 0.015;   % ~1.5 km spatial bounding box
% time_window_mins  = 45;      % +/- 45-minute temporal window
% 
% fprintf('\nStarting Matchup Extraction across %d Stations...\n', N_stations);
% 
% for i = 1:N_stations
%     st_id   = data1.St(i);
%     st_lat  = data1.Lat(i);
%     st_lon  = data1.Lon(i);
%     st_time = data1_datetime(i);
% 
%     % Spatial & Temporal Masks
%     spatial_mask  = (abs(data2_lats - st_lat) <= geo_tolerance_deg) & ...
%                     (abs(data2_lons - st_lon) <= geo_tolerance_deg);
%     time_diff_mins = abs(minutes(data2_datetime_est - st_time));
%     temporal_mask = (time_diff_mins <= time_window_mins);
% 
%     match_mask  = spatial_mask & temporal_mask;
%     num_matches = sum(match_mask);
% 
%     if num_matches > 0
%         % Calculate ensemble mean across all spectra matched within the window
%         matched_spectra = data2_rrs_array(match_mask, :);
%         matched_rrs_matrix(i, :) = mean(matched_spectra, 1, 'omitnan');
%         fprintf('Station %2d: Averaged %2d matching SoRad spectra.\n', st_id, num_matches);
%     else
%         % Nearest spatial neighbor fallback if no strict temporal match exists
%         dist_deg = sqrt((data2_lats - st_lat).^2 + (data2_lons - st_lon).^2);
%         [min_dist, nearest_row] = min(dist_deg);
%         matched_rrs_matrix(i, :) = data2_rrs_array(nearest_row, :);
%         fprintf('Station %2d: No strict match within window -> Extracted nearest spatial spectrum (Distance: %.4f deg).\n', ...
%             st_id, min_dist);
%     end
% end
% 
% %% Step 7: Export Combined Matchup CSV File
% % Base station information from Data 1
% base_info_table = data1(:, {'St', 'Lon', 'Lat', 'Date', 'TimeEST', 'TSS'});
% 
% % Matched SoRad spectra table with formatted headers
% matched_rrs_table = array2table(matched_rrs_matrix, 'VariableNames', rrs_headers_data2);
% 
% % Concatenate metadata + spectra
% final_matchup_table = [base_info_table, matched_rrs_table];
% 
% % Save output CSV
% writetable(final_matchup_table, output_csv);
% fprintf('\nSuccessfully saved matchup CSV to:\n%s\n', fullfile(pwd, output_csv));




%% Updated criteria
%% Step 1: Initialize Workspace & Set File Paths
clc; clear; close all;

% Define Working Directories & File Paths
target_dir = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end

data1_path = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/St_cruise_insitu_TSS.csv';
data2_path = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad/rrs_with_rrswl_09_2_SAV.csv'; 
output_csv = fullfile(target_dir, 'St_cruise_SoRad_Matchups_Filtered_10min.csv');

%% Step 2: Load Input Data Files
fprintf('Loading Data 1 (In-Situ Station TSS): %s\n', data1_path);
opts1 = detectImportOptions(data1_path);
opts1.VariableNamingRule = 'preserve';
data1 = readtable(data1_path, opts1);

fprintf('Loading Data 2 (SoRad Spectrometer Data): %s\n', data2_path);
opts2 = detectImportOptions(data2_path);
opts2.VariableNamingRule = 'preserve';
opts2.VariableNamesLine = 1; % Row 1 as column headers
opts2.DataLines = [2 Inf];   % Data starts on Row 2
data2 = readtable(data2_path, opts2);

%% Step 3: Standardize SoRad Wavelength Headers (340..899 nm -> Rrs_340..Rrs_899)
vars2 = data2.Properties.VariableNames;
numeric_band_idx = [];
new_vars2 = vars2;

for i = 1:length(vars2)
    var_str = char(strtrim(vars2{i}));
    var_clean = regexprep(var_str, '[^0-9a-zA-Z_]', '');
    tokens = regexp(var_clean, '(\d+)', 'tokens');
    
    if ~isempty(tokens)
        wl_num = str2double(tokens{1}{1});
        if wl_num >= 340 && wl_num <= 899
            new_vars2{i} = sprintf('Rrs_%d', round(wl_num));
            numeric_band_idx = [numeric_band_idx, i]; %#ok<AGROW>
        end
    end
end

data2.Properties.VariableNames = new_vars2;
rrs_headers_data2 = new_vars2(numeric_band_idx);

% Extract numeric array of wavelengths for QC filtering
all_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), rrs_headers_data2);

fprintf('Identified %d Rrs bands in SoRad Data 2 (%s to %s).\n', ...
    length(rrs_headers_data2), rrs_headers_data2{1}, rrs_headers_data2{end});

%% Step 4: Parse Dates & Times (UTC -> Eastern Local Time UTC-4 / EDT)
vars2_clean = data2.Properties.VariableNames;

date_col_idx = find(contains(vars2_clean, 'date', 'IgnoreCase', true), 1);
time_col_idx = find(contains(vars2_clean, 'time', 'IgnoreCase', true) | ...
                    contains(vars2_clean, 'utc', 'IgnoreCase', true), 1);

date_col_name = vars2_clean{date_col_idx};
time_col_name = vars2_clean{time_col_idx};

% Data 1 Timestamps (EST / Local Time)
data1_dt_str = strcat(cellstr(string(data1.Date)), {' '}, cellstr(string(data1.TimeEST)));
data1_datetime = datetime(data1_dt_str, 'InputFormat', 'MM/dd/yy H:mm');

% Data 2 Timestamps (Convert UTC to EDT/Eastern Local Time)
data2_dt_str = strcat(cellstr(string(data2.(date_col_name))), {' '}, cellstr(string(data2.(time_col_name))));
data2_datetime_utc = datetime(data2_dt_str, 'InputFormat', 'MM/dd/yy HH:mm:ss', 'TimeZone', 'UTC');

% Convert to Eastern Time (America/New_York)
data2_datetime_edt = datetime(data2_datetime_utc, 'TimeZone', 'America/New_York');
data2_datetime_edt.TimeZone = ''; % Remove timezone tag for clean difference math

%% Step 5: Spatio-Temporal Matchup Loop with Quality Control
N_stations = height(data1);
N_rrs_bands = length(rrs_headers_data2);
matched_rrs_matrix = NaN(N_stations, N_rrs_bands);
data2_rrs_array = table2array(data2(:, rrs_headers_data2));

% Column indices for 400 - 700 nm region QC check
qc_band_idx = find(all_wls >= 400 & all_wls <= 700);

% Clean string cell array of dates in Data 1 and Data 2
data1_date_strs = cellstr(string(data1.Date));
data2_date_strs = cellstr(string(data2.(date_col_name)));

% Matching Criteria Parameters
time_window_mins = 10; % 10-minute window

fprintf('\nStarting Matchup Extraction across %d Stations...\n', N_stations);

for i = 1:N_stations
    st_id   = data1.St(i);
    st_time = data1_datetime(i);
    st_date_str = data1_date_strs{i};
    
    % 1. Robust Date Matching via String Comparison
    date_mask = strcmp(data2_date_strs, st_date_str);
    same_date_indices = find(date_mask);
    
    if isempty(same_date_indices)
        warning('Station %2d: No SoRad records found on date %s.', st_id, st_date_str);
        continue;
    end
    
    % 2. Temporal Window Check (within 10 mins)
    time_diff_mins_all = abs(minutes(data2_datetime_edt(same_date_indices) - st_time));
    within_10_sub_idx  = find(time_diff_mins_all <= time_window_mins);
    candidate_indices  = same_date_indices(within_10_sub_idx);
    
    % 3. Quality Control Filter: Exclude negative Rrs in 400-700 nm region
    valid_candidate_indices = [];
    for k = 1:length(candidate_indices)
        row_idx = candidate_indices(k);
        spec_400_700 = data2_rrs_array(row_idx, qc_band_idx);
        
        if all(spec_400_700 >= 0)
            valid_candidate_indices = [valid_candidate_indices, row_idx]; %#ok<AGROW>
        end
    end
    
    num_valid = length(valid_candidate_indices);
    
    if num_valid > 0
        % Calculate ensemble mean across valid quality-checked spectra
        matched_spectra = data2_rrs_array(valid_candidate_indices, :);
        matched_rrs_matrix(i, :) = mean(matched_spectra, 1, 'omitnan');
        
        fprintf('Station %2d (%s %s): Averaged %2d valid SoRad spectra within 10-min window.\n', ...
            st_id, st_date_str, char(data1.TimeEST(i)), num_valid);
    else
        % Fallback: Find nearest temporal valid spectrum on the same date passing QC
        valid_date_qc_indices = [];
        for k = 1:length(same_date_indices)
            row_idx = same_date_indices(k);
            spec_400_700 = data2_rrs_array(row_idx, qc_band_idx);
            if all(spec_400_700 >= 0)
                valid_date_qc_indices = [valid_date_qc_indices, row_idx]; %#ok<AGROW>
            end
        end
        
        if ~isempty(valid_date_qc_indices)
            time_diffs_qc = abs(minutes(data2_datetime_edt(valid_date_qc_indices) - st_time));
            [min_tdiff, min_idx] = min(time_diffs_qc);
            best_row = valid_date_qc_indices(min_idx);
            
            matched_rrs_matrix(i, :) = data2_rrs_array(best_row, :);
            fprintf('Station %2d (%s %s): No 10-min match -> Used nearest valid spectrum on same date (Row %d, Time diff: %.1f mins).\n', ...
                st_id, st_date_str, char(data1.TimeEST(i)), best_row, min_tdiff);
        else
            warning('Station %2d: No valid non-negative spectra found on date %s.', st_id, st_date_str);
        end
    end
end

%% Step 6: Construct Output Table & Export CSV
base_info_table = data1(:, {'St', 'Lon', 'Lat', 'Date', 'TimeEST', 'TSS'});
matched_rrs_table = array2table(matched_rrs_matrix, 'VariableNames', rrs_headers_data2);

final_matchup_table = [base_info_table, matched_rrs_table];

% Save result to CSV
writetable(final_matchup_table, output_csv);
fprintf('\nSuccessfully saved non-NaN SoRad matchup CSV to:\n%s\n', output_csv);

%% With time stamp
%% Step 1: Initialize Workspace & Set File Paths
clc; clear; close all;

% Define Working Directories & File Paths
target_dir = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad';
if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end

data1_path = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/St_cruise_insitu_TSS.csv';
data2_path = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/Sorad/rrs_with_rrswl_09_2_SAV.csv'; 
output_csv = fullfile(target_dir, 'St_cruise_SoRad_Matchups_Filtered_10min.csv');

%% Step 2: Load Input Data Files
fprintf('Loading Data 1 (In-Situ Station TSS): %s\n', data1_path);
opts1 = detectImportOptions(data1_path);
opts1.VariableNamingRule = 'preserve';
data1 = readtable(data1_path, opts1);

fprintf('Loading Data 2 (SoRad Spectrometer Data): %s\n', data2_path);
opts2 = detectImportOptions(data2_path);
opts2.VariableNamingRule = 'preserve';
opts2.VariableNamesLine = 1; % Row 1 as column headers
opts2.DataLines = [2 Inf];   % Data starts on Row 2
data2 = readtable(data2_path, opts2);

%% Step 3: Standardize SoRad Wavelength Headers (340..899 nm -> Rrs_340..Rrs_899)
vars2 = data2.Properties.VariableNames;
numeric_band_idx = [];
new_vars2 = vars2;

for i = 1:length(vars2)
    var_str = char(strtrim(vars2{i}));
    var_clean = regexprep(var_str, '[^0-9a-zA-Z_]', '');
    tokens = regexp(var_clean, '(\d+)', 'tokens');
    
    if ~isempty(tokens)
        wl_num = str2double(tokens{1}{1});
        if wl_num >= 340 && wl_num <= 899
            new_vars2{i} = sprintf('Rrs_%d', round(wl_num));
            numeric_band_idx = [numeric_band_idx, i]; %#ok<AGROW>
        end
    end
end

data2.Properties.VariableNames = new_vars2;
rrs_headers_data2 = new_vars2(numeric_band_idx);

% Extract numeric array of wavelengths for QC filtering
all_wls = cellfun(@(x) str2double(strrep(x, 'Rrs_', '')), rrs_headers_data2);

fprintf('Identified %d Rrs bands in SoRad Data 2 (%s to %s).\n', ...
    length(rrs_headers_data2), rrs_headers_data2{1}, rrs_headers_data2{end});

%% Step 4: Parse Dates & Times (UTC -> Eastern Local Time UTC-4 / EDT)
vars2_clean = data2.Properties.VariableNames;

date_col_idx = find(contains(vars2_clean, 'date', 'IgnoreCase', true), 1);
time_col_idx = find(contains(vars2_clean, 'time', 'IgnoreCase', true) | ...
                    contains(vars2_clean, 'utc', 'IgnoreCase', true), 1);

date_col_name = vars2_clean{date_col_idx};
time_col_name = vars2_clean{time_col_idx};

% Data 1 Timestamps (EST / Local Time)
data1_dt_str = strcat(cellstr(string(data1.Date)), {' '}, cellstr(string(data1.TimeEST)));
data1_datetime = datetime(data1_dt_str, 'InputFormat', 'MM/dd/yy H:mm');

% Data 2 Timestamps (Convert UTC to EDT/Eastern Local Time)
data2_dt_str = strcat(cellstr(string(data2.(date_col_name))), {' '}, cellstr(string(data2.(time_col_name))));
data2_datetime_utc = datetime(data2_dt_str, 'InputFormat', 'MM/dd/yy HH:mm:ss', 'TimeZone', 'UTC');

% Convert to Eastern Time (America/New_York)
data2_datetime_edt = datetime(data2_datetime_utc, 'TimeZone', 'America/New_York');
data2_datetime_edt.TimeZone = ''; % Remove timezone tag for clean difference math

%% Step 5: Spatio-Temporal Matchup Loop with Quality Control & Time Capture
N_stations = height(data1);
N_rrs_bands = length(rrs_headers_data2);
matched_rrs_matrix = NaN(N_stations, N_rrs_bands);
sorad_time_est_cell = cell(N_stations, 1); % Preallocate for SoRad timestamps

data2_rrs_array = table2array(data2(:, rrs_headers_data2));

% Column indices for 400 - 700 nm region QC check
qc_band_idx = find(all_wls >= 400 & all_wls <= 700);

% Clean string cell array of dates in Data 1 and Data 2
data1_date_strs = cellstr(string(data1.Date));
data2_date_strs = cellstr(string(data2.(date_col_name)));

% Matching Criteria Parameters
time_window_mins = 10; % 10-minute window

fprintf('\nStarting Matchup Extraction across %d Stations...\n', N_stations);

for i = 1:N_stations
    st_id   = data1.St(i);
    st_time = data1_datetime(i);
    st_date_str = data1_date_strs{i};
    
    % 1. Robust Date Matching via String Comparison
    date_mask = strcmp(data2_date_strs, st_date_str);
    same_date_indices = find(date_mask);
    
    if isempty(same_date_indices)
        warning('Station %2d: No SoRad records found on date %s.', st_id, st_date_str);
        sorad_time_est_cell{i} = 'N/A';
        continue;
    end
    
    % 2. Temporal Window Check (within 10 mins)
    time_diff_mins_all = abs(minutes(data2_datetime_edt(same_date_indices) - st_time));
    within_10_sub_idx  = find(time_diff_mins_all <= time_window_mins);
    candidate_indices  = same_date_indices(within_10_sub_idx);
    
    % 3. Quality Control Filter: Exclude negative Rrs in 400-700 nm region
    valid_candidate_indices = [];
    for k = 1:length(candidate_indices)
        row_idx = candidate_indices(k);
        spec_400_700 = data2_rrs_array(row_idx, qc_band_idx);
        
        if all(spec_400_700 >= 0)
            valid_candidate_indices = [valid_candidate_indices, row_idx]; %#ok<AGROW>
        end
    end
    
    num_valid = length(valid_candidate_indices);
    
    if num_valid > 0
        % Calculate ensemble mean across valid quality-checked spectra
        matched_spectra = data2_rrs_array(valid_candidate_indices, :);
        matched_rrs_matrix(i, :) = mean(matched_spectra, 1, 'omitnan');
        
        % Record average SoRad timestamp for matched spectra
        avg_sorad_dt = mean(data2_datetime_edt(valid_candidate_indices));
        sorad_time_est_cell{i} = datestr(avg_sorad_dt, 'HH:MM:SS');
        
        fprintf('Station %2d (%s %s): Averaged %2d valid SoRad spectra (SoRad Time: %s).\n', ...
            st_id, st_date_str, char(data1.TimeEST(i)), num_valid, sorad_time_est_cell{i});
    else
        % Fallback: Find nearest temporal valid spectrum on the same date passing QC
        valid_date_qc_indices = [];
        for k = 1:length(same_date_indices)
            row_idx = same_date_indices(k);
            spec_400_700 = data2_rrs_array(row_idx, qc_band_idx);
            if all(spec_400_700 >= 0)
                valid_date_qc_indices = [valid_date_qc_indices, row_idx]; %#ok<AGROW>
            end
        end
        
        if ~isempty(valid_date_qc_indices)
            time_diffs_qc = abs(minutes(data2_datetime_edt(valid_date_qc_indices) - st_time));
            [min_tdiff, min_idx] = min(time_diffs_qc);
            best_row = valid_date_qc_indices(min_idx);
            
            matched_rrs_matrix(i, :) = data2_rrs_array(best_row, :);
            sorad_time_est_cell{i} = datestr(data2_datetime_edt(best_row), 'HH:MM:SS');
            
            fprintf('Station %2d (%s %s): No 10-min match -> Used nearest valid spectrum on same date (SoRad Time: %s, Time diff: %.1f mins).\n', ...
                st_id, st_date_str, char(data1.TimeEST(i)), sorad_time_est_cell{i}, min_tdiff);
        else
            warning('Station %2d: No valid non-negative spectra found on date %s.', st_id, st_date_str);
            sorad_time_est_cell{i} = 'N/A';
        end
    end
end

%% Step 6: Construct Output Table & Export CSV
base_info_table = data1(:, {'St', 'Lon', 'Lat', 'Date', 'TimeEST'});
base_info_table.SoRad_TimeEST = sorad_time_est_cell; % Add SoRad Time column
base_info_table.TSS = data1.TSS;                     % Attach TSS column

matched_rrs_table = array2table(matched_rrs_matrix, 'VariableNames', rrs_headers_data2);

final_matchup_table = [base_info_table, matched_rrs_table];

% Save result to CSV
writetable(final_matchup_table, output_csv);
fprintf('\nSuccessfully saved non-NaN SoRad matchup CSV to:\n%s\n', output_csv);

