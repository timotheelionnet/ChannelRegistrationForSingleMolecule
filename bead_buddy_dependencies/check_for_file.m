


function [file_found, file_path] = check_for_file(project_dir, pat)
% CHECK FOR file in input dir (does not open subdirs)
file_found = 0;
file_path = [];
% extract file list
list_dir = dir(project_dir);
% get string array of file names
f_names = string( {list_dir.name} );
folders = string( {list_dir.folder} );

% check if a config file (.ini) is present
for i = 1:numel(f_names)
    if contains(f_names(i), pat, 'IgnoreCase', 1)
        file_found = 1;
        
        file_path = fullfile(folders(i), f_names(i));
        
        
        

    end
    
end

     
end