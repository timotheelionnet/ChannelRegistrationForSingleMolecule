%% Get NN matrix
% Dependencies: compute_distMatrix_from_loc.m

% Inputs: refLoc3List = file list of reference channel loc3 paths
%         queryLoc3Lists = cell array: each cell contains a file list of
%         loc3 paths to be searched for nearest neighbors
%         nnThresh = max dist in units of loc3 file to call a NN pair

% Outputs: cell array: each cell contains a matrix of 

% If the kwArg voxSize is not specified, then analysis is done in units of
% pixels

function [nn_mat] = get_nn_matrix(refLoc3List, queryLoc3Lists, nnThresh, options)
    arguments
        refLoc3List cell
        queryLoc3Lists cell
        nnThresh double
        options.units (1,1) string = 'pix';
        options.voxSize = [1, 1, 1];
        options.nDims double = 3;
        options.roiLims double = [];
    end
    roiLims = options.roiLims;
    voxSize = options.voxSize;
    nDims = options.nDims;
    
    
    %initialize matrix that will store all matched NN pairs
    allNn_data = {}; 
    %We will count how many missing NNs for each ref spot
    all_missed = zeros(height(refLoc3List), numel(queryLoc3Lists));
    
    %Loop through reg channels
    for i = 1:numel(queryLoc3Lists)

        % ref spot countnig
        nRefSpots = 0;
        
        myQueryLoc3List = queryLoc3Lists{i}; %load the ith reg channel file list

     
        
      
        %Row counter of ith NN matrix
        ctr = 1;
        for m = 1:numel(refLoc3List) %loop through FOVs

            
    %     for m = 1:1 %test just 1
            
            missed = 0; %count ref spots without nn in mth FOV

            curRefLoc3 = load(fullfile(refLoc3List{m}),'-ascii');
            curQueryLoc3 = load(fullfile(myQueryLoc3List{m}),'-ascii');

            if isempty(curRefLoc3) || isempty(curQueryLoc3)
                continue
            end
            
            %Load spots from mth FOV for ref and reg ch and convert to nm
            %If voxSize =[1,1,1], then these functions are trivial
            refLoc3nm = convert_loc_pix_to_nm(curRefLoc3, voxSize(1:nDims));
            queryLoc3nm = convert_loc_pix_to_nm(curQueryLoc3, voxSize(1:nDims));
            
            % count number of ref spots
            nRefSpots = nRefSpots + height(curRefLoc3);

            % Apply an ROI
            if ~isempty(roiLims)

                refLoc3nm = apply_ROI(refLoc3nm, roiLims);
              

                queryLoc3nm = apply_ROI(queryLoc3nm, roiLims);
               
            end

            %get all by all dist mat for mth FOV
            myDistMat = compute_distMatrix_from_loc(refLoc3nm, queryLoc3nm, nDims); 
   
            
            %Loop through each row (ie ref spot) and find nearest neighbor
            for p = 1:height(myDistMat)
                
                %Store index of nn as nnIdx
                nnIdx = find(myDistMat(p,:) == min(myDistMat(p,:)));
                
                %Is nn within our nn threshold?
                if myDistMat(p, nnIdx) <= nnThresh
    
                    % Is spot p also the NN for spot nnIdx?
                    mutualNnIdx = find(myDistMat(:,nnIdx) == min(myDistMat(:,nnIdx)));
    
                    if mutualNnIdx ==p
                        
                        %add it to the nn dataframe: refloc3nm, regLoc3nm,
                        %dr, dx dy dz
                      
                        allNn_data{i}(ctr,:) = [refLoc3nm(p,1:3) , queryLoc3nm(nnIdx,1:3), myDistMat(p, nnIdx), queryLoc3nm(nnIdx,1:3) - refLoc3nm(p,1:3)];
                    
                    
                        %Row counter of ith nn matrix
                        ctr = ctr +1;
                    end
                else
                    %This ref spot did not have a match in the reg channel
                    missed = missed +1;
                    
                end          
            end

            % Store number of unmatched ref spots for mth FOV in query ch i
            all_missed(m,i) = missed / height(myDistMat);
        end
       
       
    end    
    %Get percentage of spots that were left unmatched per FOV per registration
    %channel
    avg_missed = (sum(all_missed, 1) / height(all_missed)) * 100;
    
    for i = 1:width(avg_missed)
        
        disp(append("Average of: ", string(avg_missed(i)), "% of ref spots had no NN in Query Ch" + string(i)));
    end
    nn_mat = allNn_data;

    disp(string(nRefSpots) + 'reference FISH spots in dataset')
    
end