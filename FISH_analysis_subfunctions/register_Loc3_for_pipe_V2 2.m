%%  APPLY REGISTRATION CORRECTION
% For pipe version V1— NOT SELF-CONTAINED. 
% Exepcts fs to be loaded and apply_channels.m to have been run

% INPUt = loc3 data with channel and FOV info encoded in filename. INI file
% should contain a loc3 file type
% Output = Loc3 table in same format as input
% Option to work from fileset or just apply to a specified table

% expects the myfits matrix to be loaded with element IJ being the model to
 
%% define channels to register and output dir

try 
    exist(regModels, "var")
catch
    regModels = reg_models;
end


% Extract list of your fish channels with data to be registered
regChannels = setdiff(allFISHChannels, refChannel);

% Make output directory
outFolder = fullfile(myParams.Output.outFolder, "res", "bead_corrected_loc");

if ~exist(outFolder, "dir")
    mkdir(outFolder);
    disp("Results will be saved to");
    disp(string(outFolder));
else
    disp("Results will be saved to existing folder");
    disp(string(outFolder));
end


%% Loop through channels and load file lists

%init cell array to store  file lists for each fish channel
allFileLists = cell(1, numel(allFISHChannels)); 

%Load loc3 file lists into cell array 
for i = 1:numel(allFISHChannels)
    allFileLists{i} = fs.getFileName({'Channel'}, {allFISHChannels(i)} ,'loc','all');

    disp(append("You have loaded ", string(height(allFileLists{i})), ...
        " Images for Channel " ,string(allFISHChannels(i))))
end    

disp( ' ');
disp(append("Channel ", string(refChannel), " is your reference channel for registration"));
disp(' ');

% postReg all file lists. This will be fore downstream analysis. 
% stores the filepaths of registered loc3 files
postReg_allFL = cell(size(allFileLists));

% add corrected loc3 file type
fs.addFileType("loc_REG", []);
disp("Added file type loc_REG to fileset")
disp(" ");

%% Save the reference channel loc3 images as is but renamed

if crop == 1
    disp('ROI applied');
    disp('X: ' + string(roi_x));
    disp('Y: ' + string(roi_y));
end

% get the indices of your reference channels in the array allFISHChannels
refIdx = find(ismember(allFISHChannels, refChannel));

% loop thru ref loc files and save as is but with new name for downstream processing
for i = 1:numel(refIdx)
    myRefFl = allFileLists{refIdx(i)};
    
    % create row in the post reg all file list table
    postReg_allFL{refIdx(i)} = cell(size(myRefFl));
    
    for j = 1:numel(myRefFl)

        myRefLoc = load(string(myRefFl(j)));

        % option to apply ROI
        if crop == 1
            
            myRefLoc = apply_ROI(myRefLoc, roi);
        end
    
        % get name of current ref ch loc file
        [f,n,x] = fileparts(myRefFl(j));
    
        % save renamed loc in pix as cleanREG but not actually registered!
        saveName = fullfile(outFolder, n + "REG" + x);
        
        % add to the postReg FL
        postReg_allFL{refIdx(i)}{j} = char(saveName);
        
        save(saveName,'myRefLoc', '-ascii');

        % find where we are in table and add reg loc3 to fl
        curFl_idx = find(ismember(fs.fList.loc, myRefFl(j)));
        fs.fList.loc_REG(curFl_idx) = {char(saveName)};
        
    end
end

%% Apply bead correction and save 


% get the indices of your reference channels in the array allFISHChannels
regIdx = find(ismember(allFISHChannels, regChannels));

% loop thru ref loc files 
for i = 1:numel(regIdx)

    myRegCh = allFISHChannels(regIdx(i));
    myRegFl = allFileLists{regIdx(i)};
    curModel = regModels{myRegCh, refChannel}; 

    % create row in the post reg all file list table
    postReg_allFL{regIdx(i)} = cell(size(myRegFl));
    
    for j = 1:numel(myRegFl)
        
        % load jth loc file for ith registration channel
        myRegLoc = load(string(myRegFl(j)));
        
        % check for empties and save
        if isempty(myRegLoc)
            % get name of current reg ch loc file
            [f,n,x] = fileparts(myRegFl(j));
        
            % save renamed loc in pix!
            saveName = fullfile(outFolder, n + "REG" + x);
            
            % add to postReg FL table
            postReg_allFL{regIdx(i)}{j} = char(saveName);
            
            save(saveName,'saveNewLoc', '-ascii');
            % find where we are in table
            curFl_idx = find(ismember(fs.fList.loc, myRegFl(j)));
            fs.fList.loc_REG(curFl_idx) = {char(saveName)};

            continue
        end
            

        myRegLocNm = convert_loc_pix_to_nm(myRegLoc, voxSize);
        spot_XY = myRegLocNm(:, 1:2);

        dx = feval(curModel{1}{1}, spot_XY);
        dy = feval(curModel{2}{1}, spot_XY);
        dz = feval(curModel{3}{1}, spot_XY);

        newLocNm = [(myRegLocNm(:,1:3) - [dx dy dz]) myRegLocNm(:,4:end) ];
        saveNewLoc = convert_loc_nm_to_pix(newLocNm, voxSize);
        
        % option to apply ROI
        if crop == 1
            
            saveNewLoc = apply_ROI(saveNewLoc, roi);
        end
    
        % get name of current reg ch loc file
        [f,n,x] = fileparts(myRegFl(j));
    
        % save renamed loc in pix!
        saveName = fullfile(outFolder, n + "REG" + x);

        
        
        % add to postReg FL table
        postReg_allFL{regIdx(i)}{j} = char(saveName);
        
        save(saveName,'saveNewLoc', '-ascii');


        % find where we are in table
        curFl_idx = find(ismember(fs.fList.loc, myRegFl(j)));
        fs.fList.loc_REG(curFl_idx) = {char(saveName)};
        
    end
end

%% Plot NN dists before and after registration

myNnThresh = 1000;
nDims = 3;

refLoc3List = fs.getFileName({'Channel'},{refChannel},'loc','all');
regChLocLists = cell(1, numel(regChannels));

for i = 1:numel(regChannels)
    
    regChLocLists{i} = fs.getFileName({'Channel'}, {regChannels(i)} ,'loc','all');        
end    

% Using the get nn matrix func Nutual NN is enforced
% inputs are ref loclist, cell array of reg loc lists, nn thresh in
% relevant units, nDims, and vox size in desired units per pixel in xyz
% Output is [refx, refy, refx], [regx, regy, regz], 3d dist, xdist,ydist, zdist
% as well as percentage of spots that did not get matched to an NN. Often
% due to different thresholdng in AIRLOCALIZE

disp('Before bead correction of FISH images');
allNn_data = get_nn_matrix(refLoc3List, regChLocLists, myNnThresh, nDims=nDims, voxSize=voxSize);
disp('');

% now lets look at post reg NN data. Nutual NN is enforced

% Using the get nn matrix func
% inputs are ref loclist, cell array of reg loc lists, nn thresh in
% relevant units, nDims, and vox size in desired units per pixel in xyz
% Output is [refx, refy, refx], [regx, regy, regz], 3d dist, xdist,ydist, zdist
% as well as percentage of spots that did not get matched to an NN. Often
% due to different thresholdng in AIRLOCALIZE

% lets PLOT

for i = 1:numel(regChannels)
    postRegChLocLists{i} = postReg_allFL{regIdx(i)};
end    

disp('After bead correction of channels: ' + string(regChannels));
allNn_data_postReg = get_nn_matrix(refLoc3List, postRegChLocLists, myNnThresh, nDims=nDims, voxSize=voxSize);
disp('');


myReg_CDF = figure;
legTxt = [];
for i = 1:numel(allNn_data_postReg)
    cdfplot(allNn_data{i}(:,7))

    legTxt = [legTxt, "Ch" + string(allFISHChannels(regIdx(i))) + "vs" + string(refChannel) + "Raw"];
    hold on
    cdfplot(allNn_data_postReg{i}(:,7))
    legTxt = [legTxt, "Ch" + string(allFISHChannels(regIdx(i))) + "vs" + string(refChannel) + "Corrected"];
end
legend(legTxt, "Location","east");

cdfFigName = 'Bead Corrected NN Analysis With Thresh=' + string(myNnThresh); 
title(cdfFigName);
xlabel('3D NN Dist (nm)');
savefig(myReg_CDF, fullfile(outFolder, cdfFigName+'.fig'));


%% 

disp('~~~~~~')
disp('Bead Correction Applied')
disp('~~~~~~')


