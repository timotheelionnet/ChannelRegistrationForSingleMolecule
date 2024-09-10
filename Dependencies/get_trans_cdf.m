function [trans_matches] =  get_trans_cdf(loci, match_thresh, mutualNN_only, nDims)

    % we will find the transallelic matching cdf
    
    % get a loc file
    cur_loc = loci;
    
    % compute lower triangle dist matrix with dimensinoality specified
    d = compute_distMatrix_from_loc(cur_loc, [] , nDims);
    
    if mutualNN_only == 1
        % get closest spot for each spot in loci
        [min_i_vals, min_i_idx] = min(d, [], 2);

        min_i_vals = min_i_vals(min_i_vals <= match_thresh);
    else
    
        idx = find(d<match_thresh);

        min_i_vals = d(idx);
        
    
    end
    
    trans_matches = min_i_vals;
    

