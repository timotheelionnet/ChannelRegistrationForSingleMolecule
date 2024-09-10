%% find FIJI files

% get the jar for ij
[ij_found, ij_jar_path] = find_jar(fullfile(FIJI_path, "jars/"), 'ij-');

disp('Found FIJI at:')
disp(ij_jar_path)

% add to path
javaaddpath(ij_jar_path);

%% find FIJI scripts path

% get scripts dir explicitly
FIJI_scripts_dir = fullfile(FIJI_path, "scripts/");
% add to path
addpath(FIJI_scripts_dir)


%% get MIJI

% find MIJ filepath in our dependencies folder

filePath = matlab.desktop.editor.getActiveFilename;

[code_dir,f,~] = fileparts(filePath);

depend_dir = fullfile(code_dir, "Dependencies/MIJI/");

[mij_found , mij_path] = check_for_file(depend_dir, 'mij.jar');


% add to path
javaaddpath(mij_path);

disp('Found MIJI at:')
disp(mij_path)





%% functions


function [file_found, file_path] = check_for_file(project_dir, pat)
% CHECK FOR filefile
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
%%

function checkFilePath(filePath)
    % checkFilePath Checks if the specified file path exists.
    %   If the file path exists, it returns true.
    %   If the file path does not exist, it returns an error message.

    if isfile(filePath)
        disp('The file exists.');
    elseif isfolder(filePath)
        disp('The folder exists.');
    else
        error('Error: The specified file path does not exist.');
    end
end

function [file_found, file_path] = find_jar(project_dir, pat)
% CHECK FOR filefile
file_found = 0;
file_path = [];
% extract file list
list_dir = dir(project_dir);
% get string array of file names
f_names = string( {list_dir.name} );
folders = string( {list_dir.folder} );

% check if a config file (.ini) is present
for i = 1:numel(f_names)
    if startsWith(f_names(i), pat, 'IgnoreCase', 1) && endsWith(f_names(i), 'jar')
        file_found = 1;
        
        file_path = fullfile(folders(i), f_names(i));
        
        
        

    end
    
end

     
end



