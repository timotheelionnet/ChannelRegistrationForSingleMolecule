%%  Channel assignment and registration function mapper

% find the fit function object in the directory. Searches all subdirs
% rootdir = fileparts(channelKeyPath); 
% first serach in sample folder
filelist_samp = dir(fullfile(sample_dir, '**', '*.*'));  %get list of files 
filelist_samp = filelist_samp(~[filelist_samp.isdir]);  %remove folders from list just in case

filelist = dir(fullfile(rootdir));  %get list of files only in parent folder
filelist = filelist(~[filelist.isdir]);  %remove folders from list just in case
% find the channel key and report to user --------------------------------
% first look in sample folder
fileFound = 0;
for i = 1:numel(filelist_samp)
    if contains(filelist_samp(i).name, "channelKey")
        channelKeyPath = fullfile(filelist_samp(i).folder, filelist_samp(i).name);
        fileFound = 1;
    end
end
%%
% if no channel key in sample folder, we look for a master channel
% key in parent folder
if fileFound ~= 1
    for i = 1:numel(filelist)
        if contains(filelist(i).name, "channelKey")
            channelKeyPath = fullfile(filelist(i).folder, filelist(i).name);
            fileFound = 1;
        end
    end
end
%%

channelKey = readtable(channelKeyPath);
regModels = cell(max(max(channelKey.uM_FISH_ch), max(channelKey.uM_Bead_ch)));

% report
if fileFound == 1
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    disp("Loaded channel key from");
    disp(channelKeyPath);
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
else
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    disp("No channel key found");
    disp("Make sure there is a 'channelKey.txt' file in the same directory tree as your cfg file");
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    return % halt script if no fits are found 
end

% find registration model functions and report to user -------------------
fileFound = 0;
% search in sample dir first
for i = 1:numel(filelist_samp)
    if contains(filelist_samp(i).name, "fits.mat")
        fitsPath = fullfile(filelist_samp(i).folder, filelist_samp(i).name); 
        fileFound = 1;
    end
end

%%
% if not see if there is one in the project dir
if fileFound ~= 1
    for i = 1:numel(filelist)
        if contains(filelist(i).name, "fits.mat")
            fitsPath = fullfile(filelist(i).folder, filelist(i).name); 
            fileFound = 1;
        end
    end
end


%%

% report to user
if fileFound == 1
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    disp("Loaded registration model from");
    disp(fullfile(fitsPath));
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
else
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    disp("No registration model object found");
    disp("Make sure there is a 'fits.mat' file in the same directory tree as your cfg file");
    disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    return % halt script if no fits are found 
end
% assemble fitfunction object filepath and load into workspace
load(fitsPath)



%% developing assignment

% go through and find the bead uM idx that corresponds 
for i = 1:numel(channelKey.uM_FISH_ch)
    if channelKey.uM_FISH_ch(i) == refChannel
        beadRefChIdx = channelKey.uM_Bead_ch(i); % this is bead ch val for ref laser line
    end
end

% loop thru all laser lines
mapError=1;
for i = 1:numel(channelKey(:,"dye"))
    
%     for ith localizable laser line
    if channelKey.isLocalizable(i) == 1
        
        % get the uM FISH idx for the current laser line
        fshIdx = channelKey.uM_FISH_ch(i);
        
        

        % get the uM bead idx that is for the same laser line
        beadIdx = channelKey.uM_Bead_ch(i);

        % get the reg model that goes from current laser line to the
        % reference laser line (using bead idxs)
        myRegModel = myFits{beadIdx, beadRefChIdx};
        
        % store the appropriate fit function in terms of the FISH indexes
        regModels{fshIdx, refChannel} = myRegModel;

        if~isempty(myRegModel)
            mapError = 0;
        end
    end
end

if mapError == 1
    error('Mapping failed: Double check your FISH reference channel and channel key')
end

% Report to user
disp(' ');
disp('Registration functions have been mapped to FISH channels');
disp('cell ij contains function mapping FISH ch i to ch j')
disp(' ');
disp(regModels);
