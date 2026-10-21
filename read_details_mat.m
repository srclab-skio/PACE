%% Read details .mat file
% Specify the name of your .mat file
filename = '/Users/masud/OneDrive02/MODIS_manuscript/StCruise26/MODIS_l2/TSS_MLP_Model_5600_Train_R2_0_75_Test.mat';

% Display details of the variables inside the file (Name, Size, Bytes, Class)
fprintf('Details for %s:\n', filename);
fileDetails = whos('-file', filename);

% Convert the struct array to a table for a cleaner display
detailsTable = struct2table(fileDetails);

% Display the table 
% (Note: Some older MATLAB versions might require displaying the struct directly)
disp(detailsTable(:, {'name', 'size', 'bytes', 'class'}));

% --- Optional: If you want to actually load the data ---
% Load the file contents into a struct variable
data = load(filename);
disp(data);

