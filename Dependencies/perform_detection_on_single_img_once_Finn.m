%% Finn Clark 6/18/2024
% This function allows you to call AIRLOCALIZE 2 for a given image file and
% allows you set parameters in-line. Only the first three args are
% required. Output is loc file and par4 filepath

function [loc, par4] = perform_detection_on_single_img_once_Finn(imgpath, interactivePSF, interactiveThresh, psfXY, psfZ, threshType, threshVal)

% store the template string inside this litttle for loop
for d =1
    my_par4 =     ["[Files]";
        "dataFileName=img_file_path";
        "saveDirName=save_dir_path";
        "fileProcessingMode=singleFile";
        "inclusionString=";
        "recursive=0";
        "";
        "[detectionSettings]";
        "psfSigma=psfXY, psfZ";
        "threshUnits=SD";
        "threshLevel=4";
        "maxSpots=20000";
        "fitMethod=3DMaskFull";
        "";
        "[userInput]";
        "setPsfInteractively= set_PSF_interactively";
        "setThreshInteractively=set_thresh_interactively";
        "movieFrameUsedToSetParameters=1";
        "outputSpotsImage=0";
        "autoSpotSize=1";
        "spotSize=2, 2";
        "autoSpotIntensity=1";
        "spotIntensity=1000";
        "";
        "[advanced]";
        "minDistBetweenSpots=2";
        "filterLo=2.15";
        "filterHi=1";
        "fittedRegionSize=3";
        "tol=0.01";
        "bgRegionThickness=1";
        "psfType=integratedGaussian";
        "maxIterations=100";
        "bgCorrectionMode=localPlane";
        "useLegacyFormats=0";
        "";
        ""];
end

my_img_filepath = imgpath;

% PSF interactivitiy and defaults
if interactivePSF == 1 
    % interactivitiy 
    user_set_PSF = string(1);

        
    
else
    
    user_set_PSF = string(0);
    psfXY = string(1);
    psfZ = string(2);

    % psf vals
    my_par4 = strrep(my_par4, 'psfXY', string(psfXY));
    my_par4 = strrep(my_par4, 'psfZ', string(psfZ));

end

% Threshold interactiveity
if interactiveThresh == 1 
    % interactivitiy 
    user_set_thresh = string(1);
else
    user_set_thresh = string(0);

    if isempty(threshType) || isempty(threshVal)
        error('You must specify thresh type and value: types are ''SD'' or ''Abs'' ');
              
    end

    my_par4 = strrep(my_par4, "threshLevel=4", strcat("threshLevel=",string(threshVal)) );
end

% load generic params file





% savedir
[my_save_dir,save_prefix,~] = fileparts(my_img_filepath);








%%
%____________Modify the par4 file____________

% modify image file path 
my_par4 = strrep(my_par4, 'img_file_path', my_img_filepath);
%%
% modify save dir (default is same as image data
my_par4 = strrep(my_par4, 'save_dir_path', my_save_dir);


% user interactivitiy
my_par4 = strrep(my_par4, 'set_PSF_interactively', user_set_PSF);
my_par4 = strrep(my_par4, 'set_thresh_interactively', user_set_thresh);




save_name = strcat(save_prefix, '.par4');
% save(fullfile(my_save_dir, save_name), "my_par4", '-ascii')

% save the new par4
new_par4_filepath = char(fullfile(my_save_dir, save_name));
writelines(my_par4, new_par4_filepath);

loc = AIRLOCALIZE(new_par4_filepath);
par4 = new_par4_filepath;
end








