

function [fits_found, fits_path] = check_for_fits(project_dir)
% check for fits file
[fits_found, fits_path] = check_for_file(project_dir, 'fits.mat');


% report to user 
if fits_found == 1
    disp('~~~~~')
    disp('fits.mat found at:')
    disp(fits_path)
else
    disp('~~~~~')
    disp('no fits.mat file present in:');
    disp(project_dir)

end
end
