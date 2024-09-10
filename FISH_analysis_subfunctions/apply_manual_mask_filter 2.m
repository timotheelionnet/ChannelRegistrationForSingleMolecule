
% name of file type to be filtered
file2filter = 'loc3_REG_sorted';


% load manual mask file

% must have allFISHCHannels loaded

user_maskT = readtable("G:\Finn\20231009_PPARG_FluorSwap_trial2\day14_rep2_1\adipo_manual_mask.csv");

labels = user_maskT.Label; % we will have to use some regex stuff or pattern matching on this

mask_x = user_maskT.X; 

mask_y = user_maskT.Y; % remember FIJI records Y from the top down... WHY! we will correct this below

% base_str is the common string preceding the unique FOV number in your
% mask table

base_str = 'rep2_1_';

%%

% loop thru data by FOV

my_FOVs = unique(fs.fList.FOV);


for i = 1:numel(my_FOVs)

    cur_FOV = my_FOVs(i);

    % loop thru each channel for current FOV 
    for j = 1:numel(allFISHChannels)
%     for j =1:1

        cur_ch = allFISHChannels(j);

        % get required files
        cur_filepath = fs.getFileName({'FOV', 'Channel'}, {cur_FOV, cur_ch}, file2filter);

        cur_loc = load(cur_filepath{1});

        cur_mask_path = fs.getFileName({'FOV', 'Channel'}, {cur_FOV, cur_ch}, 'mask');

        cur_mask_img = readTifStackWithImRead(cur_mask_path{1});

        % get the x and y from our manual filter
        % the are the indices of the rows relevant to cur loc file
        keep_idx = find( contains(labels, strcat(base_str, string(cur_FOV)) ));
        
        % units of pixels! THE X Y FROM FIJI MUST BE SWITCHED BRO
        keep_coords = [user_maskT.Y(keep_idx), user_maskT.X(keep_idx)];
        

        % column index number where cellIDs are stored
        cellID_idx = width(cur_loc) + 1;
        


        % sort the loci
        
        sorted_loc = sort_loci_into_cell_masks(cur_loc, cur_mask_img);

        % firure out cellIDs of the cells we clicked on
        
        
        cells2keep = sort_loci_into_cell_masks(keep_coords, cur_mask_img); %% the cellID will be in col 3
        

        if ~isempty(cells2keep)
            % remove loc that are not in desired cells
        
            rmv_idx = ~(ismember(sorted_loc(:, cellID_idx), cells2keep(:, 3)));
    
            cur_loc(rmv_idx, :) = [];
        else
            % there are no cells to keep
            cur_loc(:,:) = [];

        end

        [f,n,x] = fileparts(cur_filepath{1});
        saveName = fullfile(outFolder, strcat(n, x));
        save(saveName,'cur_loc', '-ascii');

        disp('saved')
        saveName
        cur_loc


        % find where we are in table and add reg loc3 to fl
        % first need to find col idx of the file type we are sorting
        curFileType_idx = find(ismember(fs.fList.Properties.VariableNames, file2sort));
        % then we find what row we are in
        tempFl = fs.fList(:,curFileType_idx);
        curFl_row_idx = find(ismember(tempFl{:,1}, curLocFl{j}));
        
        addFileFl_idx =  find(ismember(fs.fList.Properties.VariableNames, string(file2sort)));
        % ADD function to SAVE
        fs.fList(curFl_row_idx, addFileFl_idx) = {char(saveName)};

        





       

        

    end

end








%%

% load 