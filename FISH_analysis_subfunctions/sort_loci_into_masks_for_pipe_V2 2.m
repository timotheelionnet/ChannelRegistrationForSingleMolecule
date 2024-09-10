%% Sorting spots into cell masks and saving as new loc file

% for modularity we should load our fileset and reload allFISH channels
% we also should have the name of the file type we are sorting as an input


% add column to file type for the sorted spot files
fs.addFileType(file2sort + "_sorted", []);
disp("~~~~~~~~~~~~~~~~~~~~~~");
disp("Added file type " + file2sort + "_sorted" + " to FS");
disp(" ");

% check for res subfolder
outFolder = fullfile(myParams.Output.outFolder, "res", "sorted_spots");
if ~exist(outFolder, "dir")
    mkdir(outFolder);
    disp("Results will be saved to");
    disp(string(outFolder));
else
    disp("Results will be saved to existing folder");
    disp(string(outFolder));
    disp(" ");
    disp("~~~~~~~~~~~~~~~~~~~~~~");
end

%% loop through all FISH img files and sort based on associated mask files

% QC raw spots per cell
spot_qc = cell(size(allFISHChannels));
spot_qc2 = cell(size(allFISHChannels));

for i = 1:numel(allFISHChannels)
    
    % loop thru each fish CH
    curCh = allFISHChannels(i);
    disp('sorting spots in channel ' + string(curCh));

    if apply_spot_thresh == 1
        disp(string(spotThresh) + 'spot thresh applied');
    end

    % get all mask files
    curMaskFl = fs.getFileName({'Channel'}, {curCh}, {'mask'}, {'all'});

    % get all loc3 files to sort
    curLocFl = fs.getFileName({'Channel'}, {curCh}, {file2sort}, {'all'});

    % err check
    if numel(curLocFl) ~= numel(curMaskFl)
        disp("You do not have a mask file for every loc file.")
        disp("Check that your mask files are named according to FISH Channel and FOV values ")
        return
    end
    
    % go through each loc file and sort into mask
    % save spots per cell per FOV before and after filtering
    spot_qc_FOV = [];
    spot_qc_FOV2 = [];
    for j = 1:numel(curLocFl)
        

        myLoc = load(curLocFl{j}, '-ascii');

        if isempty(myLoc)
            [f,n,x] = fileparts(curLocFl{j});
            saveName = fullfile(outFolder, n + "_sorted" +x);
            save(saveName,'mySortedLoc', '-ascii');
    
    
            % find where we are in table and add reg loc3 to fl
            % first need to find col idx of the file type we are sorting
            curFileType_idx = find(ismember(fs.fList.Properties.VariableNames, file2sort));
            % then we find what row we are in
            tempFl = fs.fList(:,curFileType_idx);
            curFl_row_idx = find(ismember(tempFl{:,1}, curLocFl{j}));
            
            addFileFl_idx =  find(ismember(fs.fList.Properties.VariableNames, string(file2sort+"_sorted")));
            % ADD function to SAVE
            fs.fList(curFl_row_idx, addFileFl_idx) = {char(saveName)};
            continue
        end

        myMask = readTifStackWithImRead(curMaskFl{j});
        cellID_indx= width(myLoc) +1;
        % adds column with cell ID. ID = 0 means in the background!
        mySortedLoc = sort_loci_into_cell_masks(myLoc, myMask);
        
        % figure
        % imagesc(myMask)
        % hold on
        % scatter(myLoc(:,2), myLoc(:,1), '*', 'red')
        
        % option to remove bkrnd spots
        if rmv_bkrnd_spots == 1
            bkrndSpot_idx = find(ismember(mySortedLoc(:, width(mySortedLoc)), 0));

            mySortedLoc(bkrndSpot_idx,:) = [];
%             disp(string(numel(bkrndSpot_idx)) + " Background spots with cellID = 0 removed");
        end
        
        % get spots per cell for QC
        [spotsPerCell, myIDs] = get_spots_per_cell(mySortedLoc, cellID_indx);
        % get spots per cell for each FOV and add to list for this channel
        spot_qc_FOV = [spot_qc_FOV; spotsPerCell];

        % option to gate by number of cells
        if apply_spot_thresh == 1
            
            cellID_indx=7;
    
            % spotsPerCell is column with numel = unique cells in loc, vals = spots per unique cellID
            %myIDs is cellID value in loc file
            [spotsPerCell, myIDs] = get_spots_per_cell(mySortedLoc, cellID_indx);


        
        
    
            % Remove cells w more than spotThresh spots
            %find index of cells with too many spots, ie fals positives (fp)
            fpCellIdx = find(spotsPerCell < spotThresh(1) | spotsPerCell > spotThresh(2));
            % get the cellID values for cells
            fpID = myIDs(fpCellIdx);
            % remove spots from loc1 corresponding to those cells

%             disp(mySortedLoc);
            rmvIdx = find(ismember(mySortedLoc(:,cellID_indx), fpID));
            mySortedLoc(rmvIdx, :) = [];
%             disp(mySortedLoc);
            [spotsPerCell, myIDs] = get_spots_per_cell(mySortedLoc, cellID_indx);
        end
        % we will also save the histogram of spots per cell after filtering
        spot_qc_FOV2 = [spot_qc_FOV2; spotsPerCell];

        % save the sorted loc3 files

        [f,n,x] = fileparts(curLocFl{j});
        saveName = fullfile(outFolder, n + "_sorted" +x);
        save(saveName,'mySortedLoc', '-ascii');


        % find where we are in table and add reg loc3 to fl
        % first need to find col idx of the file type we are sorting
        curFileType_idx = find(ismember(fs.fList.Properties.VariableNames, file2sort));
        % then we find what row we are in
        tempFl = fs.fList(:,curFileType_idx);
        curFl_row_idx = find(ismember(tempFl{:,1}, curLocFl{j}));
        
        addFileFl_idx =  find(ismember(fs.fList.Properties.VariableNames, string(file2sort+"_sorted")));
        % ADD function to SAVE
        fs.fList(curFl_row_idx, addFileFl_idx) = {char(saveName)};


    end
    
    % the ith cell contains spots per cell for all j FOVs
    spot_qc{i} = spot_qc_FOV;
    spot_qc2{i} = spot_qc_FOV2;




end

% histogram for raw spots per cell no filtering done for these
for i = 1:numel(spot_qc)
    raw_spot_fig = figure;
    curSpots = spot_qc{i};
    
    histogram(curSpots);
    xlabel('n spots per cell')
    ylabel('n')
    title('Ch' + string(allFISHChannels(i))+ 'Raw Spots Per Cell')
    savefig(raw_spot_fig, fullfile(outFolder, 'Ch' + string(allFISHChannels(i))+'raw_spot_per_cell_hist'));

    thresh_spot_fig = figure;
    

    curSpots2 = spot_qc2{i};
    
    histogram(curSpots2);
    xlabel('n spots per cell')
    ylabel('n')
    title('Ch' + string(allFISHChannels(i))+ ' Filtered Spots Per Cell')
    
    
    savefig(thresh_spot_fig, fullfile(outFolder, 'Ch' + string(allFISHChannels(i)) + 'filtered_spot_per_cell_hist'));
%     map = tab10; 
%     figure
%     histogram(curSpots, 'FaceAlpha', .1)
%     hold on
%     histogram(curSpots2, 'FaceAlpha', .2)
     
%     box off
%     axis tight
%     
%     legend boxoff
end

