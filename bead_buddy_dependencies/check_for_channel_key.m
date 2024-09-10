

function [key_found, key_path] = check_for_channel_key(project_dir)
% check for channel key file
[key_found, key_path] = check_for_file(project_dir, 'key');


% report to user 
if key_found == 1
    disp('~~~~~')
    disp('channel key found at:')
    disp(key_path)
else
    disp('~~~~~')
    disp('no channel key present in:');
    disp(project_dir)

end
end