
function bead_dir = find_bead_dir(directory)
    % findSubdirectoryContainingBead Searches for a subdirectory containing 'bead' (ignoring case).
    %   The function returns the path to the subdirectory if found.
    %   If no such subdirectory is found, it returns an empty array.

    % Get a list of all files and folders in the specified directory
    items = dir(directory);
    
    % Filter out the items to keep only directories
    isDir = [items.isdir];
    subDirs = items(isDir);
    
    % Initialize the output variable
    bead_dir = '';
    
    % Loop through each subdirectory to check for the presence of 'bead' in the name
    for i = 1:length(subDirs)
        subDirName = subDirs(i).name;
        % Skip the current and parent directory entries
        if strcmp(subDirName, '.') || strcmp(subDirName, '..')
            continue;
        end
        % Check if 'bead' is present in the subdirectory name, ignoring case
        if contains(lower(subDirName), 'bead')
            bead_dir = fullfile(directory, subDirName);
            return; % Exit the function as soon as the first match is found
        end
    end
    
    % If no matching subdirectory is found, return an empty array
    if isempty(bead_dir)
        disp('No subdirectory containing "bead" was found.');
    end
end