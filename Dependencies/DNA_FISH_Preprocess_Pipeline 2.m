%% DNA FISH PRE-PROCESSING PIPELINE (run one chunk at a time)
% requires
% -MIJI install
% -FIJI with biovoxxel plug in
% -Add Ons: Computer Vision, Deep Learning, GPU Coder, Medical Imaging,
% Medical Image Interface for Cellpose

% ~~~~~~~Before starting, make sure that:~~~~~~~
% --You have added the pipeline directory to your path with all subfolders
% --all FOVs per condition are in the same folder 
% ----'beads' is not in FISH data folder name

% --your bead stacks are in a folder with the word 'beads' in it


% %% USER INPUT: Project directory
% % give the path to the project folder containing bead and/or FISH image folders
% project_dir = 'H:\FC_DNA_FISH_Pipeline_Install_20240618\20240618_Pipeline_Export_AL2\demo';

%% USER INPUT: Pre-process Bead images

% nChannels = 4; % number of channels in ur bead hyperstacks
% doFlatField =1; % apply flat field correction to your iamges
% 
% % ref_channel = 2; % the channel you will map your data to. Use uManager channel value in your beads, eg 'C3-beads.loc3'
% 
% voxSize = [98, 98, 250]; % x,y,z of your voxel in nm
% 
% % whether to run localization
% runAirlocalize = 1;
%     % set to 1 to localize and analyze bead images
%     % set to 0 if loc files already present


%% Find MIJI and FIJI

% FIJI_path = 'C:\Users\Lionnet Lab\Documents\Fiji.app';

% code_dir = fileparts(matlab.desktop.editor.getActiveFilename);
% cd(code_dir)
% 
% % Generate a path string including all subfolders
% pathWithSubfolders = genpath(code_dir);
% 
% % Add the path to MATLAB's search path
% addpath(pathWithSubfolders);
% 
% % Display a message indicating that the path has been added
% disp(['Added the following directory to the MATLAB path: ' code_dir]);
% 
% 
% %~~~~RUN~~~~%
% MIJ_installer
% 
% checkFilePathForSpaces(project_dir);
% checkFilePathForSpaces(code_dir);

%% CHANNEL KEY
% check for channel_key

[key_found, key_path] = check_for_channel_key(project_dir);

% generate channel key template test GUI
if ~key_found
    % key_path = generate_key_tab(project_dir, 5);
    generate_channel_key_gui_V2
    
    disp('please fill out the table template via GUI')
    [key_found, key_path] = check_for_channel_key(project_dir);
    
end


% load key
key_tab = readtable(key_path);

% check for errors in key
verify_channel_key(key_tab)

% interpret channel key and write a cfg file
[~,parent] = fileparts(project_dir);
ref_ch = key_tab.uM_FISH_ch( key_tab.isReference == 1);
bead_ref_ch = key_tab.uM_Bead_ch( key_tab.isReference == 1);
bead_ch_vals = key_tab.uM_Bead_ch( key_tab.uM_Bead_ch ~=0);


key_tab 


%% BEADS
% if we are using bead correction
if use_bead_buddy == 1
    %~~~~RUN~~~~~~%
    [bead_dir, bead_img_dir] = organize_bead_imgs(project_dir);
end

%% WRITE CFG BEADS

% if we are using bead correction
if use_bead_buddy == 1
    bead_paths = get_tif_list(bead_img_dir);
    
    % make split channel dir
    split_ch_dir = fullfile(bead_dir, 'split_channels');
    if ~exist (split_ch_dir, 'dir')
        mkdir(split_ch_dir)
        
        disp('~~~~~')
        disp('bead imgs will be split and saved to:')
        disp(split_ch_dir)
    end
    
    % check for config file
    % [cfg_found, cfg_path]  = check_for_cfg(project_dir);
    
    % write a cfg file 
    cfg_path = bead_cfg_writer_for_pipe(code_dir, split_ch_dir, parent, bead_ch_vals, bead_ref_ch);
end

%% MIJI BEADS

% if we are using bead correction
if use_bead_buddy == 1
    %~~~RUN~~~%
    miji_process_bead_hyperstacks(bead_paths, split_ch_dir);
end

%% AIRLOCALIZE BEADS

% if we are using bead correction
if use_bead_buddy == 1
    %~~~RUN~~~%
    MODULAR_Bead_Analysis_Pipeline_V6 
end
%% Organize bead correction functions

% parentPath = 'D:\Finn\20240229_Staggered_PPARG_DNAFISH_Anlaysis\subFOlder'; % user defined multi or single dataset

if use_bead_buddy == 1
    %~~~RUN~~~%
    organize_bead_buddy_functions_for_pipe
end
%% END OF BEAD BUDDY GRAFT


%%  USER INPUT: Pre-process FISH images

% mask_chan = 1; % specify uManager channel index in FISH images to be used for segmentation
% do_maxZ = 1; % save max z projections of the mask_chans?
% make_cell_masks = 1; % make cell masks with cellpose? Requires do_maxZ =1
% make_cfg = 1; % save a config file for this sample? (Required for analysis)
% 
% process_images = 1; %do you want to actually process the images? Set 0 to only organize files / or make masks / write cfg
% doFlatField =1; % apply flat field correction (I always use this)
% 
% % For cell segmentation
% avg_cell_diameter= 150; % avg nuclear diameter in pixels
%     % 70 pix for B-cells at 100x, 150 pix for MSCs at 60x Nikon


    % You can use your own cell masks but they must be named according to
    % the formula: 'conditionName_FOV.tif_ZMAX_MASK.tif
        % yes I realize this kinda sucks
%%
%~~~~RUN~~~~~~%
 % NEED TO MAKE IT SO IF PROCSES IMAGES IS OFF WE DONT NEED THIS
FISH_channels = key_tab.uM_FISH_ch(key_tab.isLocalizable == 1);

% get rid of the stuff that tries to check for cfg and image org
Integrated_MIJI_preProcess_v2

disp('')
disp('~~~~~')
disp('Pre-process pipeline complete')
disp('Proceed to analysis')
disp('~~~~~')

    