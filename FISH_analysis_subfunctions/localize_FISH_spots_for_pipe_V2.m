% Localize spots for PIPE
% 5/22/2023

% Check for cfg file
try
    disp(config_fname);
catch
    disp('Remember to specify a path to cfg file');
end

% Convert to char
config_fname = char(config_fname);

% create file list object
fs = fileSet;

fs.buildFileSetFromConfig(config_fname);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','Date','FOV'},convertNonNumeralStrings);

% read parameters from config file relative to analysis into a parameters object
myParams = readParamsFromIniFile(config_fname);

outFolder = fullfile(myParams.Output.outFolder, "res", "localization");
if ~ exist(outFolder,'dir')
    mkdir(outFolder);
end

voxSize = [myParams.Settings.voxSize_dxy,...
            myParams.Settings.voxSize_dxy,...
            myParams.Settings.voxSize_dz];

%% Airlocalize

%Add a column to your fileSet object to hold localization data
fs.addFileType('loc',[]);

%Extract ist of your channels with data to be localized
myChannels = myParams.channelDescription.fishChannelList; 

%Initialize cell array that will contain a column of image filepaths for each
%channel. 
myFileLists = cell(1, numel(myChannels)); 


%Loop through your channels and extract the filenames for all images
%in that channel. Display file list size 

for i=1:numel(myParams.channelDescription.fishChannelList) %index based on how many fish channels
    myFileLists(i) = {fs.getFileName({'Channel'},{myChannels(i)},'fishImg','all')};   

    disp(append("You have loaded ", string(height(myFileLists{i})), ...
    " FISH Images for Channel " ,string(myChannels(i))))
end
% Print ref channel

disp(append("Channel ", string(refChannel), " is your reference channel for spot matching"));
%% GUI Airlocalize
trainAIRLOCALIZE = 1; % True = Train and then run AIRLOCALIZE on iamge data. False = skip the manual training/thresholding step (you already have loc3 files saved).
if trainAIRLOCALIZE == 1
    interactive_thresh = 1;
    interactive_PSF = 0; % change if for some reason this matters
    
    % -----------------------------GUI for getting PSF and thresh

    
    %I want to go through and manually get fit params using thefirst image of
    %each channel.
     
    locpar = cell(1, numel(allFISHChannels));
    for i = 1:numel(allFISHChannels)
        disp('Processing');
        disp(myFileLists{i}{1});
        % [~, par4_path] = perform_detection_on_single_img_once_Finn(myFileLists{i}{1}, 0, trainAIRLOCALIZE,1,2); %Manual training on first image for ith channel

        [~, par4_path] = perform_detection_on_single_img_once_Finn(myFileLists{i}{1}, interactive_PSF, interactive_thresh, 1,2); %Manual training on first image for ith channel
        
        locpar{i} = par4_path; %Store first par4 file for ith channel
       
        
    end
end
%% If you have your fit params already in your workspace, run from here

if trainAIRLOCALIZE ==1 || runAIRLOCALIZE == 1
    %Now I have trained AIRLOCALIZE on all channels and stored the fit params
    %Time to loop through all of my images and analyze those using previosuly
    %stored fit params
    imageCount= append("ch:", string(numel(myChannels)), " file:", string(numel(myFileLists{1})));
    
    
    for i = 1:numel(myChannels)
        for j=1:numel(myFileLists{i})
            
            
            disp(append("processing file ch:", string(i), " file:", string(j), " out of ", imageCount));
    
            
            %Using ith channel parameters, replace data with jth image in
            %channel i
            % locpar{i}.data = readTifStackWithImRead(myFileLists{i}{j});
    
            % run airlocalize without GUI for the rest of the files
            loc = parlooper(myFileLists{i}{j}, locpar{i});
            
            disp(append("Processed", string(myFileLists{i}{j})));
            
            %For each image we want to store loc4 data and its file name in the
            %correct row of our loc4 fileSet column. loc4 is overwritten each
            %iteration
            % loc4 = res{i}.final_pix;
    
            % build the loc filename to enter in the fileSet list
            [d,f,~] = fileparts(myFileLists{i}{j});
            locFileName = fullfile(d,[f,'.loc4']);
            
            % find the conditions of the current file
            fl = fs.getFileName({},{},'fishImg','all');
            idx = find(ismember(fl,myFileLists{i}{j}));
            
          
            % curCond = {};
            % 
            % for k = 1:numel(fs.conditions)
            %     curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
            % end
            % %Save the images loc4 filename in the fileSet
            % fs.setFileName(fs.conditions,curCond,'loc4',locFileName);
            
            %Now save the actual file
            save(locFileName,'loc','-ascii');
            
    
        end
    end

else
    disp('Not running AIRLOCALIZE');
    disp('');
end

%% Reload fs to get loc files in fs flist
fs = fileSet;

fs.buildFileSetFromConfig(config_fname);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','Date','FOV'},convertNonNumeralStrings);

%% loop through FISH Channels and plot dist hist
nDims = numel(voxSize);

% compute a histogram of randomly localized spots for histogram normalization
% find image size
curChannel = myParams.channelDescription.fishChannelList(1);
curImgFileName = fs.getFileName({'Channel'},{curChannel},'fishImg','first');
curImgFileName = curImgFileName{1};
imSize = getImageSizeFromFile(curImgFileName); % get image size in pixels




%% All by all normalized dist hists by channel
% convert to nm
% voxSize = [myParams.Settings.voxSize_dxy,...
%             myParams.Settings.voxSize_dxy,...
%             myParams.Settings.voxSize_dz];
% imSize_nm = convert_loc_pix_to_nm(imSize,voxSize);
% 
% [nRand,randBinsCenters] = generate_pairDist_hist_randSpots(...
%     myParams.Settings.nRandSpots,imSize_nm,...
%     myParams.Settings.distMatrixHistogram_binSize);
% 
% % generate output folder for histogram figures
% outFolderHist = fullfile(outFolder,'pairDistHist');
% if ~ exist(outFolderHist,'dir')
%     mkdir(outFolderHist);
% end
% 
% for i=1:numel(myParams.channelDescription.fishChannelList)
%     curChannel = myParams.channelDescription.fishChannelList(i);
%     curFileList = fs.getFileName({'Channel'},{curChannel},'loc3','all');
%     for j=1:numel(curFileList) 
%         
%         disp('Processing file:')
%         [~,myF,~] = fileparts(curFileList{j});
%         disp(myF);
%             
%         % load loc file
%         loc = load(curFileList{j},'-ascii');
% 
%         if isempty(loc)
%             continue
%         end
%         
%         % convert to nm (assumes 3D)
%         loc3nm = convert_loc_pix_to_nm(loc,voxSize);
%         
%         % compute distance matrix
%         [dr, ~] = compute_distMatrix_from_loc(loc3nm,[],nDims);
%         
%         % compute pair distance histogram
%         [n,x] = hist(dr( ~isnan(dr) ), randBinsCenters);
%         
%         % normalize
%         fh = figure('name','Normalized pair distances between spots');
%         plot(x,n./nRand);
%         xlabel('Distance in nm');
%         ylabel('Number of spot pairs');
%         
%         % save histogram
%         [~,f,~] = fileparts(curFileList{j});
%         savefig(fh,fullfile(outFolderHist,[f,'_pairDistHist.fig']));
%         close(fh);
%     end
% end

%% loop through FISH Channels and clean up double detections
% load image info

if clean_double_detections == 1

    disp('Removing Double Detections');
    
    voxSize = [myParams.Settings.voxSize_dxy,...
                myParams.Settings.voxSize_dxy,...
                myParams.Settings.voxSize_dz];

    nDims = numel(voxSize);
    
    curChannel = myParams.channelDescription.fishChannelList(1);
    curImgFileName = fs.getFileName({'Channel'},{curChannel},'fishImg','first');
    curImgFileName = curImgFileName{1};
    imSize = getImageSizeFromFile(curImgFileName); % get image size in pixels
    imSize_nm = convert_loc_pix_to_nm(imSize,voxSize);
    
    % fs.addFileType('cleanLoc3',[]);
    
    % Make a subdir to store reconstructed gauss spots
    gaussSpotDir= fullfile(outFolder, 'reconsructed_gauss_spots');
    if ~exist(gaussSpotDir,'dir')
        mkdir(gaussSpotDir);
    end
    
    % Loop through FISH channels and clean up double detections
    for i=1:numel(myParams.channelDescription.fishChannelList)
        curChannel = myParams.channelDescription.fishChannelList(i);
        curFileList = fs.getFileName({'Channel'}, {curChannel},'loc','all');
        
        if isempty(curFileList)
            disp('No loc3 files in your fileSet. Make sure they exist and then reload your fileset');
        end
    
        for j=1:numel(curFileList) 
            disp('Processing file:')
            [~,myF,~] = fileparts(curFileList{j});
            disp(myF);
            
            % load loc file
            loc = load(curFileList{j},'-ascii');
    
            if isempty(loc)
                continue
            end
            
            % convert to nm (assumes 3D)
            loc3nm = convert_loc_pix_to_nm(loc,voxSize);
            
            % clean up spots: remove double detectin by setting min dist
            % between localizations. 
            cleanSpots = clean_up_double_detection(loc3nm, myParams.Settings.distThreshold, nDims);
            
            % output cleaned spot list
            [d,f,~] = fileparts(curFileList{j});
            
            %make sub dir for clean loc3 so they dont pollute FS building
            cleanDir = fullfile(d,'clean_loc');
            if ~exist(cleanDir, 'dir')
                mkdir(cleanDir)
            end
            
            % save the cleaned spots
            locFileName = fullfile(cleanDir,[f,'_clean.loc3']);
            pos = convert_loc_pix_to_nm(cleanSpots,1./voxSize);
            save(locFileName,'pos','-ascii');
            
            % if you want to generate gauss spots
            if gen_gauss_spots == true
                imSize2 = [imSize(2),imSize(1),imSize(3)];
                spotSize = [myParams.Settings.spotSizeInImg,myParams.Settings.spotSizeInImg,myParams.Settings.spotSizeInImg];
                s = generate_gauss_spots_in_stack(pos,spotSize,imSize2);
                save_as_tiff(s,fullfile(gaussSpotDir,[f,'_spots.tiff']));
            end
            
    
            % find the conditions in the current file
            fl = fs.getFileName({},{},'loc','all');
            idx = find( ismember(fl,curFileList{j}) );
            curCond = {};
            for k = 1:numel(fs.conditions)
                curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
            end
            fs.setFileName(fs.conditions,curCond,'loc',locFileName);
            
        end
    end
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp("Spot localization complete");
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~');