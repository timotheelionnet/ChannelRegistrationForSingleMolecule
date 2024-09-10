function [n,x] = generate_pairDist_hist_randSpots(nRandomSpots,imSize,binSize)
    % generates nRandSpots spots randomly distributed in space in a 3D
    % volume with 
    % max X = imSize(1)
    % max Y = imSize(2)
    % max Z = imSize(3)
    % outputs the histogram of pair distances (n,x) where x are the centers
    % of the bins;
    
    randSpots = rand(nRandomSpots,3);
    
    randSpots(:,1) = randSpots(:,1)*imSize(1);
    randSpots(:,2) = randSpots(:,2)*imSize(2);
    randSpots(:,3) = randSpots(:,3)*imSize(3);
    
    
    [dr, ~] = compute_distMatrix_from_loc(randSpots,[],3);
    
    
    
    [n,x] = hist(dr( ~isnan(dr) ), 0:binSize:max(dr(:)));
    
end