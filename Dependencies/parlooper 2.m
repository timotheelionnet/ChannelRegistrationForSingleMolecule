%% Finn Clark 12/4/2023
% this function uses AIRLOCALIZE 2 to process an image specified by imfilepath
% using the localization parameters written in the .par4 file specified by
% refpar4

function loc = parlooper(imgfilepath, refpar4)

% extract file name info
[my_save_dir,save_prefix,~] = fileparts(imgfilepath);


% read the reference par4 file as lines 
my_par4 = readlines(refpar4);

% loop thru and replace the image file path 
for i = 1:numel(my_par4)

    curLine = my_par4{i};

    if contains(curLine, 'dataFileName=')

        curLine = strcat('dataFileName=', char(imgfilepath));
        my_par4{i} = curLine;

        my_par4 = regexprep(my_par4,'setPsfInteractively=1', 'setPsfInteractively=0');
        my_par4 = regexprep(my_par4,'setThreshInteractively=1', 'setThreshInteractively=0'); 
        break % found it, break loop
    end

end

% save new par4 and run airlocalize
save_name = strcat(save_prefix, '.par4');
% save(fullfile(my_save_dir, save_name), "my_par4", '-ascii')

% save the new par4
save_name = strcat(save_prefix, '.par4');

new_par4_filepath = char(fullfile(my_save_dir, save_name));
writelines(my_par4, new_par4_filepath);

% run AIRLOCALIZE

loc = AIRLOCALIZE(new_par4_filepath);

end