% Clear command window
clc; clear all;

% Assume Rrs_SoRad_PACE_OLCI is already loaded
% If not, load it: load('/Users/masud/OneDriveUGA/CruiseRVSavannah/Subpixel/your_file.mat');

% Define station numbers (excluding St13)
stations = [10, 11, 12, 14, 15, 16];

% Check if required fields exist and have correct dimensions
if ~isfield(Rrs_SoRad_PACE_OLCI, 'Rrs_Sorad') || size(Rrs_SoRad_PACE_OLCI.Rrs_Sorad, 2) ~= 560
    error('Rrs_Sorad is missing or has incorrect dimensions (expected [N x 560]).');
end
if ~isfield(Rrs_SoRad_PACE_OLCI, 'Rrs_PACE') || size(Rrs_SoRad_PACE_OLCI.Rrs_PACE, 2) ~= 172
    error('Rrs_PACE is missing or has incorrect dimensions (expected [N x 172]).');
end
if ~isfield(Rrs_SoRad_PACE_OLCI, 'Rrs_OLCI') || size(Rrs_SoRad_PACE_OLCI.Rrs_OLCI, 2) ~= 16
    error('Rrs_OLCI is missing or has incorrect dimensions (expected [N x 16]).');
end

% Initialize Stations structure
Rrs_SoRad_PACE_OLCI.Stations = struct();

% Number of rows per station (modify as needed)
numRows = 5; % Assuming 1 row per station; change to 5 for multiple rows or provide external data

% Add station fields
for i = 1:length(stations)
    stationNum = stations(i);
    fieldName = sprintf('St%d', stationNum); % Create field name (St10, St11, etc.)
    
    % Initialize station structure
    Rrs_SoRad_PACE_OLCI.Stations.(fieldName) = struct();
    
    % Assign data for Sorad, PACE, and OLCI
    if i <= size(Rrs_SoRad_PACE_OLCI.Rrs_Sorad, 1)
        % Sorad: [numRows x 560]
        Rrs_SoRad_PACE_OLCI.Stations.(fieldName).Sorad = ...
            Rrs_SoRad_PACE_OLCI.Rrs_Sorad(i:min(i+numRows-1, size(Rrs_SoRad_PACE_OLCI.Rrs_Sorad, 1)), :);
    else
        Rrs_SoRad_PACE_OLCI.Stations.(fieldName).Sorad = NaN(numRows, 560);
        warning('Row %d not available for Sorad in %s; initialized with NaN.', i, fieldName);
    end
    
    if i <= size(Rrs_SoRad_PACE_OLCI.Rrs_PACE, 1)
        % PACE: [numRows x 172]
        Rrs_SoRad_PACE_OLCI.Stations.(fieldName).PACE = ...
            Rrs_SoRad_PACE_OLCI.Rrs_PACE(i:min(i+numRows-1, size(Rrs_SoRad_PACE_OLCI.Rrs_PACE, 1)), :);
    else
        Rrs_SoRad_PACE_OLCI.Stations.(fieldName).PACE = NaN(numRows, 172);
        warning('Row %d not available for PACE in %s; initialized with NaN.', i, fieldName);
    end
    
    if i <= size(Rrs_SoRad_PACE_OLCI.Rrs_OLCI, 1)
        % OLCI: [numRows x 16]
        Rrs_SoRad_PACE_OLCI.Stations.(fieldName).OLCI = ...
            Rrs_SoRad_PACE_OLCI.Rrs_OLCI(i:min(i+numRows-1, size(Rrs_SoRad_PACE_OLCI.Rrs_OLCI, 1)), :);
    else
        Rrs_SoRad_PACE_OLCI.Stations.(fieldName).OLCI = NaN(numRows, 16);
        warning('Row %d not available for OLCI in %s; initialized with NaN.', i, fieldName);
    end
end

% Verify the updated structure
disp('Updated Rrs_SoRad_PACE_OLCI structure:');
fieldNames = fieldnames(Rrs_SoRad_PACE_OLCI);
for i = 1:length(fieldNames)
    varName = fieldNames{i};
    if strcmp(varName, 'Stations')
        stationFields = fieldnames(Rrs_SoRad_PACE_OLCI.Stations);
        for j = 1:length(stationFields)
            stName = stationFields{j};
            fprintf('Stations.%s.Sorad: %s\n', stName, mat2str(size(Rrs_SoRad_PACE_OLCI.Stations.(stName).Sorad)));
            fprintf('Stations.%s.PACE: %s\n', stName, mat2str(size(Rrs_SoRad_PACE_OLCI.Stations.(stName).PACE)));
            fprintf('Stations.%s.OLCI: %s\n', stName, mat2str(size(Rrs_SoRad_PACE_OLCI.Stations.(stName).OLCI)));
        end
    else
        fprintf('%s: %s\n', varName, mat2str(size(Rrs_SoRad_PACE_OLCI.(varName))));
    end
end

% Optional: Plot Rrs for a sample station (e.g., St10)
if isfield(Rrs_SoRad_PACE_OLCI.Stations, 'St10')
    figure;
    subplot(3, 1, 1);
    plot(Rrs_SoRad_PACE_OLCI.wl_Sorad, Rrs_SoRad_PACE_OLCI.Stations.St10.Sorad(1, :), '-o');
    xlabel('Wavelength (nm)'); ylabel('Rrs (sr^-1)'); title('St10: SoRad Rrs');
    grid on;
    
    subplot(3, 1, 2);
    plot(Rrs_SoRad_PACE_OLCI.wl_PACE, Rrs_SoRad_PACE_OLCI.Stations.St10.PACE(1, :), '-o');
    xlabel('Wavelength (nm)'); ylabel('Rrs (sr^-1)'); title('St10: PACE Rrs');
    grid on;
    
    subplot(3, 1, 3);
    plot(Rrs_SoRad_PACE_OLCI.wl_OLCI, Rrs_SoRad_PACE_OLCI.Stations.St10.OLCI(1, :), '-o');
    xlabel('Wavelength (nm)'); ylabel('Rrs (sr^-1)'); title('St10: OLCI Rrs');
    grid on;
end

% Save the updated structure
save('updated_Rrs_SoRad_PACE_OLCI.mat', 'Rrs_SoRad_PACE_OLCI');