function cleanSpots = clean_up_double_detection(loc3nm, threshold, nDims)
    %get dist matrix for all spots in a loc file
    [dr, ~] = compute_distMatrix_from_loc(loc3nm,[],nDims);
    toKeep = [];
    %loop through rows of dist mat and find cols where dist is less than
    %threshold, ie will be called as a double detection
    for i=1:size(dr,1)
        %store indices of spots within the threshold
        spotsWithinThreshold = find( dr(i,:) < threshold );
        % if there are no double detections for spot i
        if isempty(spotsWithinThreshold)
            % and if we have not already stored this spot 
            if ~ismember(i,toKeep)
                % then add the low file row index i of the spot to our
                % matrix
                toKeep = [toKeep;i];
            end
        else
            % else there are multiple spots within the thresh
            % load the intensities for each spot within the threshold
            intensities = loc3nm(spotsWithinThreshold, nDims+1);
            % get the index j of the brightest spot 
            [~,j] = max(intensities);
            % if we have not already stored this spot in toKeep
            if ~ismember(j,toKeep)
                %add its index to our matrix
                toKeep = [toKeep;j];
            end
        end
    end
    %output is the input loc3nm file filtered to only include the spots
    %whose indices we stored in toKeep
    cleanSpots = loc3nm(toKeep,:);

end