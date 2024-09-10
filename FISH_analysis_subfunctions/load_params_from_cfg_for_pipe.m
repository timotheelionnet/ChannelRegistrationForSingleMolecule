% Load general info from the cfg files we generated in the preprocessing steps

config_fname = char(config_fname);
myParams = readParamsFromIniFile(config_fname);

sample_dir = myParams.Output.outFolder;
rootdir = fileparts(config_fname);

fs = fileSet;

% populate object properties from ini config file
fs.buildFileSetFromConfig(config_fname);
fs.toNumeric({'Channel','threshVal','FOV'},1);

%Read the parameters for your imaging session/data
myParams = readParamsFromIniFile(config_fname);

% load voxel dims from ini file
voxSize = [myParams.Settings.voxSize_dxy, ...
    myParams.Settings.voxSize_dxy,...
    myParams.Settings.voxSize_dz];

%Store all FISH channels (ie the C(#) in a variable
allFISHChannels = myParams.channelDescription.fishChannelList;

% Report fList size
numF = numel(unique(fs.fList.FOV));
numCh = numel(unique(fs.fList.Channel));
disp("~~~~~~Loaded Fileset Containing~~~~~~");
disp(string(numF) +" FOVs");
disp(string(numCh) + " Channels");
disp(" ");