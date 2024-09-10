function loc = sort_loci_into_cell_masks(loc,refImg)
    %take 2D or 3D loc3 and sort spots into 2D cell masks
    
    %check if loc is empty
    if ~isempty(loc)
        % into ROIs as NaN in a newly created column

        %Check if spot is outside of the image
        outOfBoundsIdx = ...
            loc(:,1) <= 0 | loc(:,1) > size(refImg,1) | ...
            loc(:,2) <= 0 | loc(:,2) > size(refImg,2);
        loc(outOfBoundsIdx,end+1) = NaN; %adds a column and sets row val to NaN if out of image bounds
        
        % Converts cell mask img, ie refImg, to linear pixel array and
        % finds the pixel our loc spots fall within
        linIdx = sub2ind(size(refImg),...
            ceil(loc(~outOfBoundsIdx,1)),...
            ceil(loc(~outOfBoundsIdx,2)));
                
        % The pixel value of the mask image is evaluated at the 
        % pixel index of our spot. If it is in a mask, then the pixel value will
        % correspond to a unique value in that FOV, if outside of a mask,
        % it will be set to 0
        
        loc(~outOfBoundsIdx,end) = refImg(linIdx);
        
    else
        disp('loc3 file is empty!')
    end
end