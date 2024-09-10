function [dr, dxyz]= compute_distMatrix_from_loc(loc1,loc2,nDims)
    
% if loc1 and loc2 are non empty matrices, returns the full distance matrix
% where dr_ij = euclidian distance between loc1(i,:) and loc2(j,:)

% if loc1 is non-empty and loc2 is empty, returns the lower triangle matrix
% where dr_ij = euclidian distance between loc1(i,:) and loc1(j,:)
% dr = compute_dist_from_loc(loc1,[],nDims)

%set nDims to 2 or 3 for 2D or 3D distances.

%%
% test if computing a distance on a single set of spots or two distinct
% sets
if isempty(loc2)
    loc2 = loc1;
    isSelfDistance = 1;
else
    isSelfDistance = 0;
end


% compute loc1 x loc2 distances
for di = 1:nDims
    dxyz{di} = ...
        repmat(loc1(:,di),1,size(loc2,1))...
        - repmat(loc2(:,di)',size(loc1,1),1);

    %remove redundant and self pairs if correlating channel
    %with itself
    if isSelfDistance
        sd = size(dxyz{di},1);
        dxyz{di}(tril(true(sd))) = NaN;
    end

    if di == 1
        dr = dxyz{di}.^2;
    else
        dr = dr + dxyz{di}.^2;
    end 
end 


dr = sqrt(dr);

end