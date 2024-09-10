%% FISH Pre process script
% IMPORTANT: your FISH condition folders CANNOT have the word 'bead' in them!!!


% % add path to Fiji matlab scripts
% addpath("C:\Users\Lionnet Lab\Documents\Fiji.app\scripts");
% % download mij.jar and move to matlab jar folder
% % make sure you have ij.jar 
% % add java paths
% javaaddpath("C:\Program Files\MATLAB\R2023b\java\jar\ij-1.54f.jar");
% javaaddpath("C:\Program Files\MATLAB\R2023b\java\jar\mij.jar");
% run miji!
Miji();
MIJ.help %tells u what commands u can do
% which system used?

% add path to dependencies
%%
% project_dir = 'D:\Finn\20240229_Staggered_PPARG_DNAFISH_Anlaysis\subFOlder'; % user defined multi or single dataset
% mask_chan = 1; % specify uManager channel index to be used for masking (these will not be flatfield corrected)
% do_maxZ = 1; % make max z projections of the mask_chans
% make_cell_masks = 1;
% make_cfg = 1;
% 
% process_images = 1;
% doFlatField_data =1; 
% voxSize = [98, 98, 250]; % x,y,z of your voxel in nm
% 
% avg_cell_diameter= 150; 
% % avg_cell_diameter= 70; b cells
% cp_thresh = 0;


%% Search for fish img directories

all_FISH_dirs = {}; % this will contain raw images to split and process
sampleDirs = {};  % thsi sill contain names of unique samples


subDirCtr = 0;
inImgDir = 0; % goes to 1 when we found the sub dir with images
% find folders in project_dir called FISHimgs

all_files = dir(fullfile(project_dir, '**/*.*'));

all_files = struct2table(all_files);

% get sub table with only images
img_idx = contains(all_files.name, 'tif');
img_table = all_files(img_idx, :);

% get unique directory names that contain images
dataDirs = unique(img_table.folder);


%%
% remove beads and any cases wehre split_channels has already been made
bead_rmv_idx = [];
for i = 1:numel(dataDirs)

    if contains(string(dataDirs{i}), 'bead', 'IgnoreCase',true) || contains(string(dataDirs{i}), 'split_channels', 'IgnoreCase',true) || contains(string(dataDirs{i}), 'max_z', 'IgnoreCase',true)
        disp('removing unwanted dir');
        
     
        bead_rmv_idx = [bead_rmv_idx; i]; % get cond dir to report to user
    end


end

dataDirs(bead_rmv_idx) = [];



disp('');
disp(string(numel(dataDirs)) + ' Unique image folders were found:');
disp(string(dataDirs));
disp('~~~~~~~~~');
disp('');

%%
% check if we already made FISH_imgs dir and cfgs
% skip_file_org = 0;
% skip_cfg = 0; 

% we will check if FISH img dirs have already been made
% bad_strings = ["FISH_imgs", "split_channels", "max_z", "res", "cell_masks"];

% rmv_idx = [];
for i = 1:numel(dataDirs)


    % check for bad strings to remove
    cur_dataDir = dataDirs{i};

    % has_bad_string = any(arrayfun(@(substr) contains(cur_sampleDir, substr, 'IgnoreCase', true), bad_strings));
    % 
    % if has_bad_string
    %     rmv_idx = [rmv_idx; i];
    % end


    % if the data dir is already in the FISH_imgs dir
    if contains(cur_dataDir, "FISH_imgs")
        disp('FISH_imgs dir already present');
        % skip_file_org = 1;

        all_FISH_dirs{i,1} = cur_dataDir; % save the FISH imgs dir
        sampleDirs{i} = fileparts(cur_dataDir); % get sample dir to report to user
    else
        
        % store the dir to be made
        all_FISH_dirs{i,1} = fullfile(cur_dataDir, 'FISH_imgs');
        sampleDirs{i} = cur_dataDir;

    end


end
% sampleDirs(rmv_idx) = [];


disp('');
disp(string(numel(sampleDirs)) + ' Unique samples were found:');
disp(string(sampleDirs'));
disp('~~~~~~~~~');
disp('');

%%


% create a FISH img dir inside each unique sample dir if not present

for i = 1:numel(all_FISH_dirs)

    disp('');   
    disp('~~~~~~~~~');

    disp('Processing files in dir: ');
    disp(string(all_FISH_dirs{i}));
    disp('');
    
    % this may or may not be the same as the FISH imgs dir
    cur_dataDir = dataDirs{i};
    
    % curFISHdir may or may not have been made yet
    curFISHdir = all_FISH_dirs{i}; 
  

    if ~exist(curFISHdir, 'dir')
        mkdir(curFISHdir);
    end

    % if contains(sampleDirs{i}, 'FISH_imgs')
    %     curFISHdir = sampleDirs{i};
    % else
    % 
    %     curFISHdir = fullfile(sampleDirs{i}, 'FISH_imgs');
    % end
    
    % make the fish imgs dir if it does not exist

    
    % add FISH dir to list if not already present

    % all_FISH_dirs = vertcat(all_FISH_dirs, curFISHdir);
    

    
    
    
    % get images for current sample

    my_FOVpaths = img_table(strcmp(img_table.folder, cur_dataDir), 'name');


   
    % define root name for images
    [~,root_name] = fileparts(sampleDirs{i});
    root_name = string(root_name);
    % disp('Renaming files as: ' + string(root_name));
    % disp('Moving files to: ' + string(curFISHdir));
    disp('');   
    disp(string(numel(my_FOVpaths.name)) + ' FOVs in this dataset');

    for j = 1:numel(my_FOVpaths.name)
        % give it a name with FOV index at the end
        new_name = root_name + '_' + string(j) + '.tif';

        old_path = fullfile(dataDirs{i}, my_FOVpaths.name{j});
        new_path = fullfile(curFISHdir, new_name);
        
        % move all indexed FOVs to the FISH img folder
        if old_path ~= new_path
            movefile(old_path, new_path);
        end
    end

%     loop thru all images in current sample folder
    
% skip organizatin if already done
    % if ~contains(curFISHdir, 'FISH_imgs')
    %     disp('Creating FISH img dirs and renaming files');
    %     for j = 1:numel(my_FOVpaths.name)
    %         % give it a name with FOV index at the end
    %         new_name = root_name + '_' + string(j) + '.tif';
    % 
    %         old_path = fullfile(sampleDirs{i}, my_FOVpaths.name{j});
    %         new_path = fullfile(curFISHdir, new_name);
    % 
    %         % move all indexed FOVs to the FISH img folder
    %         movefile(old_path, new_path);
    %     end
    % end

    

end

disp('~~~~~');
disp('File Organization Complete');
disp('~~~~~');
%% loop thru sample folders and split stacks + do z projections if desired



disp('');
disp('~~~~~~~~~');
disp('Splitting multi-channel images');
disp('');



% loop thru FISH img dirs one condition at a time
for i = 1:numel(all_FISH_dirs)

    [curSampleDir,~] = fileparts(all_FISH_dirs{i});

    % make split channel dir
    split_ch_dir = fullfile(curSampleDir, 'split_channels');
    if ~exist (split_ch_dir, 'dir')
        mkdir(split_ch_dir)
    end

    % make maxZ dir
    if do_maxZ ==1
        max_z_dir = fullfile(curSampleDir, 'max_z');
        if ~exist (max_z_dir, 'dir')
            mkdir(max_z_dir);
        end

    end


    % get all FOVS in cur fish dir FOVs
    cur_FOVs = dir(fullfile(all_FISH_dirs{i}, '*.tif'));

    cur_FOVs = struct2table(cur_FOVs);
    if process_images == 1
        % loop thru and split
        for j = 1:numel(cur_FOVs(:,"name"))
    
            if numel(cur_FOVs(:,"name")) == 1
    
                my_FOV = cur_FOVs.name;
            
                my_FOVdir = cur_FOVs.folder;
        
                my_FOVpath = fullfile(string(my_FOVdir) , string(my_FOV));
            
            else
           
                my_FOV = cur_FOVs.name(j);
                
                my_FOVdir = cur_FOVs.folder(j);
        
                my_FOVpath = fullfile(string(my_FOVdir) , string(my_FOV));
            end
            %MIJI open my_FOVpath 
            
            disp('');
            
            disp('Processing')
            disp(my_FOVpath);
            
            % close all open windows
            MIJ.closeAllWindows();
            
            % open the image in FIJI
            MIJ.run('Open...', strcat('path=', my_FOVpath));
            
            % get the title
		    curTitle = char(MIJ.getCurrentTitle());
    
            MIJ.run('')
            
            % split channels
            MIJ.run('Split Channels');
            
            splitImg_list = MIJ.getListImages;
            
            nCh = numel(splitImg_list);
            
            cl = 1:nCh;
    
            FISH_channels = cl(~ismember(cl, mask_chan));
    
            for m=1:nCh
                
                if m == mask_chan
                    
                    MIJ.selectWindow(strcat('C', string(m), '-', curTitle));
                    
                    saveTitle = char(MIJ.getCurrentTitle);
                    % save the split ch
                    MIJ.run('Save', strcat('Tiff..., path=[', fullfile(split_ch_dir, saveTitle), ']'));
    
                    if do_maxZ ==1
                                    
                        % do the max z and save that
                        
                        MIJ.selectWindow(strcat('C', string(m), '-', curTitle));
                        MIJ.run('Z Project...', 'projection=[Max Intensity]');
        %                 MIJ.selectWindow(strcat('C', string(m), '-', curTitle));
                        
                        saveTitle = strcat(saveTitle, '_ZMAX.tif');
                        MIJ.run('Save', strcat('Tiff..., path=[', fullfile(max_z_dir, saveTitle), ']'));
                        
                        disp('');
                        disp('max_z saved to');
                        disp('');
                        disp(fullfile(max_z_dir, saveTitle));
                    end
                    
                else
                    % not the mask channel, so we just split, flatfield
                    % correct, and save
                    
                    MIJ.selectWindow(strcat('C', string(m), '-', curTitle));
                    
                    if doFlatField_data ==1 
                        MIJ.run('Pseudo flat field correction', 'blurring=50 hide stack');
                    end
                    
                    saveTitle = char(MIJ.getCurrentTitle);
                    MIJ.run('Save', strcat('Tiff..., path=[', fullfile(split_ch_dir, saveTitle), ']'));
                   
                end
    
    
            end    
        end % end of ch loop

        
    end

    
    
end

% close all open windows
MIJ.closeAllWindows()
MIJ.exit

%% cfg files
for i = 1:numel(all_FISH_dirs)

    [curSampleDir,~] = fileparts(all_FISH_dirs{i});
    % write a cfg file
    if make_cfg == 1

        %   make a cfg file
        [~,root_name] = fileparts(curSampleDir);
        FISH_cfg_writer_for_pipe(curSampleDir, root_name, FISH_channels, voxSize)
        disp('cfg file written to ' + string(fileparts(all_FISH_dirs{i})));
    end

    % make a channel key for user to manipulate
    % generate_channelKey(project_dir, nCh)

end

%% cell maasks
disp('~~~~~')
disp('MAKING CELL MASKS')
disp('~~~~~')
for i = 1:numel(all_FISH_dirs)

    [curSampleDir,~] = fileparts(all_FISH_dirs{i});
    % generate cell masks
    if make_cell_masks == 1

        max_z_dir = fullfile(curSampleDir, 'max_z');

        mask_dir = fullfile(curSampleDir, 'cell_masks');
        if ~exist (mask_dir, 'dir')
            mkdir(mask_dir);
        end

        save_dir = mask_dir;

        my_fl = getFileListFromDir(max_z_dir, 0, 0);

        
        
        for p = 1:numel(my_fl)
        
            cur_zmax_path = my_fl{p};
        
            generate_cellpose_mask(cur_zmax_path, avg_cell_diameter, cp_thresh, save_dir);
        end
    end


end



%%
% MIJ.closeAllWindows();
% MIJ.exit();

    
disp('~~~~~')
disp('FISH Img processing complete')
disp('~~~~~')

    

