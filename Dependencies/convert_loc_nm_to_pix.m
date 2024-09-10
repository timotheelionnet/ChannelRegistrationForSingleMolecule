function locOut = convert_loc_nm_to_pix(locIn,voxSize)
    
    % voxSize = [dx dy] for 2D images is the size of a voxel in nm
    % voxSize = [dx dy dz] for 3D images is the size of a voxel in nm
    
    if (numel(voxSize) ~= 2) && (numel(voxSize) ~= 3)
        disp('Cannot convert loc positions to nm; voxel Size has the wrong number of elements');
        locOut = [];
        return
    end
    
    if isempty(locIn)
        disp('Loc3 file is empty!')
        return
    end
    
    locOut = locIn;
    locOut(:,1) = locIn(:,1)/voxSize(1);
    locOut(:,2) = locIn(:,2)/voxSize(2);
    if numel(voxSize) == 3
        locOut(:,3) = locIn(:,3)/voxSize(3);
    end
end