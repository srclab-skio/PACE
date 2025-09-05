%% Data subset index-based.
% Read both CSV files
data1 = readtable('GLORIA_Rrs_02Sept2025_all_Vars_350_900_TSS.csv');
data2 = readtable('Filtered_Subset_global_GLORIA_LAB_close_tss5.csv');

% Extract GLORIA_ID values from both datasets
gloria_ids_data1 = data1.GID;
gloria_ids_data2 = data2.GLORIA_ID;

% Find the indices of data1 where GLORIA_ID matches those in data2
[is_member, idx] = ismember(gloria_ids_data1, gloria_ids_data2);

% Create the subset containing only rows from data1 that have matching GLORIA_ID in data2
subset_data = data1(is_member, :);

% Display information about the result
fprintf('Original data1 size: %d rows x %d columns\n', size(data1));
fprintf('Filtered subset size: %d rows x %d columns\n', size(subset_data));
fprintf('Number of matching GLORIA_ID values: %d\n', sum(is_member));

% Optional: Save the subset to a new CSV file
writetable(subset_data, 'Filtered_Subset_GLORIA_Rrs_all_Vars_350_900_TSS.csv');
fprintf('Subset saved to: Filtered_Subset_GLORIA_Rrs.csv\n');

% Display first few rows of the subset (optional)
disp('First few rows of the subset:');
disp(head(subset_data));