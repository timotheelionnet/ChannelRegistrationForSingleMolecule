

function [cfg_found, cfg_path] = check_for_cfg(project_dir)
% check for .ini file
[cfg_found, cfg_path] = check_for_file(project_dir, 'ini');


% report to user 
if cfg_found == 1
    disp('~~~~~')
    disp('cfg file found at:')
    disp(cfg_path)
else
    disp('~~~~~')
    disp('no cfg file present in:');
    disp(project_dir)

end
end
