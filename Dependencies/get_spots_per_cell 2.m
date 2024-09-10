

function [spotsPerCell, cellID] = get_spots_per_cell(loc3, cellID_indx)

    %get all cell IDs for this match file
    myCellIDs = loc3(:, cellID_indx);
     % extract unique IDs, c == ID vals, ic = index in c of cellID value from myCellIDs
    [c, ia, ic] = unique(myCellIDs);

    % counts is how many times each unique value corresponding to indices ic is 
    % found in our loc3 file
    counts = accumarray(ic, 1);
    % store cellID value and corresponding spots per cell
    val_counts = [c, counts];

    %exclude any spots with cellID = 0
    val_counts = val_counts(val_counts(:,1)~=0, :);

    spotsPerCell = val_counts(:,2);
    cellID = val_counts(:,1);
end