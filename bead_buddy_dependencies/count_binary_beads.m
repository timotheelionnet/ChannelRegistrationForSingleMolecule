%% binarizing bead images

% input is a max z of beads
% im_path = "G:\Finn\20230901Devin_ecDNA_Rosettes\rosette_WT\max_z\C1-rosette_WT_1.tifMAX_Z.tif";
im_path = "H:\Finn_DriveH\20240226_PPARG_Staggered_Timecourse\beads\session2_beads\split_channels\C1-session2_beads_2.tif";

[I, cmap] = imread(im_path);

% crop 
% I = I(700:1000, 700:1000);

imagesc(I);

%%get a file set
config_fname = "G:\Finn\20231009_bead_binarization\bead_binary.ini";
config_fname = char(config_fname);

fs = fileSet;

% populate object properties from ini config file
fs.buildFileSetFromConfig(config_fname);
fs.toNumeric({'Channel','FOV'},1);


%Read the parameters for your imaging session/data
myParams = readParamsFromIniFile(config_fname);


%%

Igray = mat2gray(I);

imshow(Igray);


%%
%convert into grayscale image. 
k = Igray;
%calculate threshold using Otsu's method. 
level=graythresh(k); 
  
%convert into binary image using level. 
k1=imbinarize(k,level); 
  
%display the binarized image. 
imtool(k1); 


[L, num] = bwlabel(k1);

num


%% airlocalize loop until we get right number of spots

myChannels = myParams.channelDescription.fishChannelList; 

for i=1:numel(myParams.channelDescription.fishChannelList) %index based on how many fish channels
    myFileLists(i) = {fs.getFileName({'Channel'},{myChannels(i)},'fishImg','all')};   

    disp(append("You have loaded ", string(height(myFileLists{i})), ...
    " FISH Images for Channel " ,string(myChannels(i))))
end

test_file = myFileLists{i}{1};
%%
AIRLOCALIZE("G:\Finn\20231009_bead_binarization\test_file.ini")
% AL needs an ini file in the workspace in order to take image path as an
% argument
% perform_detection_on_single_image_once_Dipankar2


