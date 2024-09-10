
function [cdf_corrected_table] = cdf_purifier(my_cdf_table, trans_matches)

%     my_cdf_table = readtable("G:\Finn\20231009_PPARG_FluorSwap_trial2\day14_rep2_1\res\matched_spots\cdfs12_Sp3to6_Ch3to4_inCellOnly_noOffsectCorxn.txt");
%     my_cdf_table = readtable("G:\Finn\20231009_PPARG_FluorSwap_trial2\day0_rep1_2\res\matched_spots\cdfs12_Sp3to6_Ch3to4_inCellOnly_noOffsectCorxn.txt");
%     
    % lets make the cdf purifier into a function
    
    
    obs_cdf_x = my_cdf_table.cdfR_x;
    obs_cdf_f = my_cdf_table.cdfR_f;
    
    
    
    
    % plot cdf before analysis
%     figure
%     plot(obs_cdf_x, obs_cdf_f)
%     trans_cdf = cdfplot(trans_matches{2});
    trans_cdf = cdfplot(trans_matches);
    
    
    % need the trans data too
%     figure
    % trans_cdf = [trans_matches{2}; trans_matches{3}];
%     trans_cdf = cdfplot(trans_matches);
    
    
    trans_cdf_x = trans_cdf.XData;
    trans_cdf_f = trans_cdf.YData;
    
    
    % trans_cdf1 = cdfplot(trans_matches{2});
    % trans_cdf2 = cdfplot(trans_matches{3});
    % 
    % 
    % trans_cdf1_x = trans_cdf.XData;
    % trans_cdf1_f = trans_cdf.YData;
    % 
    % 
    % trans_cdf2_x = trans_cdf.XData;
    % trans_cdf2_f = trans_cdf.YData;
    
    %%
    % x = [trans_matches{2}; trans_matches{3}];
    % 
    % figure
    % cdfplot(x)
    
    %% interpolate both
    
    % add point 0,0 to trans data for interpolation
    trans_cdf_x = [0, trans_cdf_x];
    trans_cdf_f = [0, trans_cdf_f];
    
    % i have to sample the two curves along the same interval first, then do
    % the FIT!!!
    
    % the domain for the fit
    xq = (1:5000)';
    
    
    % remove nans
    xv = obs_cdf_x( ~( isnan(obs_cdf_x)  | isinf( obs_cdf_x ) ));
    yv = obs_cdf_f( ~( isnan(obs_cdf_x)  | isinf( obs_cdf_x ) ));

    xv = obs_cdf_x( ~( isnan(obs_cdf_f)  | isinf( obs_cdf_f ) ));
    yv = obs_cdf_f( ~( isnan(obs_cdf_f)  | isinf( obs_cdf_f ) ));
    
    % get the unique vals of y + indices then filter x, then vice versa
    [yv, ia, ic] = unique(yv);
    xv = xv(ia);
    
    [xv, ia, ic] = unique(xv);
    yv = yv(ia);
    
%     figure 
%     plot(xv, yv)
    obs_interp = interp1(xv, yv, xq);
    %%
    % --------------------------repeat with the trans data
    
    % remove nans and infs
    xv = trans_cdf_x( ~( isnan(trans_cdf_x)  | isinf( trans_cdf_x ) ));
    yv = trans_cdf_f( ~( isnan(trans_cdf_x)  | isinf( trans_cdf_x ) ));

    xv = trans_cdf_x( ~( isnan(trans_cdf_f)  | isinf( trans_cdf_f ) ));
    yv = trans_cdf_f( ~( isnan(trans_cdf_f)  | isinf( trans_cdf_f ) ));
    
    % get the unique vals of y + indices then filter x, then vice versa
    [yv, ia, ic] = unique(yv);
    xv = xv(ia);
    
    [xv, ia, ic] = unique(xv);
    yv = yv(ia);
    
    % get the interpolated values (this is an arraay not a functino output);
    % trans_interp = interp1(xv, yv, xq);
    
    % try tims code for poly fitting
    [c] = fit_3rd_degree_poly_no_constant2(xv, yv,5000);

    ytim = c(1)*xv.^3 + c(2)*xv.^2 + c(3)*xv; 

    figure
    plot(xv ,yv);
    hold on
    plot(xv, ytim);
    title('tim fit')
      
    trans_interp =  c(1)*xq.^3 + c(2)*xq.^2 + c(3)*xq; 


%     p_trans_interp = polyfit(xv, yv, 3);
%     trans_interp = polyval(p_trans_interp, xq);

    % do not allow negative values
    trans_interp(trans_interp <0) = 0;
    
    % smooth the trans
    
    % trans_interp = smooth(trans_interp, 1000);
    
    
    figure
    plot(xq, trans_interp);
    hold on
    plot(xq, obs_interp);
    hold on
    plot(trans_cdf_x, trans_cdf_f);
    title('Interplating  the CDFS')
    legend({'trans interp', 'obs interp', 'raw trans'});
    
    %% fit 
    
    % domain over which functinos are fit
    fit_domain = [750, 4500];
    
    obs2fit = obs_interp( fit_domain(1):fit_domain(2) );
    %%
    trans2fit = trans_interp( fit_domain(1):fit_domain(2) );

%     figure
%     plot(fit_domain(1):fit_domain(2), obs2fit)
%     hold on
%     plot(fit_domain(1):fit_domain(2), trans2fit)

    nan_rmv = isnan(obs2fit) | isnan(trans2fit);

    obs2fit(nan_rmv) = [];
    trans2fit(nan_rmv) = [];
    
    % analytical solution to the chi square minimization problem
    num = sum( trans2fit.^2 - (obs2fit .* trans2fit)  - trans2fit +  obs2fit );
    denom = sum( (trans2fit - 1).^2 );
    
    f= num/denom;
%     'f is'
%     num
%     denom
    disp(f)
    %%
    x = 1:5000;
    
    % find where we have nans in the trans itnerp
    
    nan_idx = isnan(trans_interp(x));
    
    % replace with 0s 
    trans_interp_crxn = trans_interp(x);
    
    % if nan, dont alter the observed
    trans_interp_crxn(isnan(trans_interp_crxn)) = 0;
    
    % generate the corrected cis matching data
    cis_corrected = (obs_interp(x) -  ((1-f) * trans_interp_crxn) ) / f;
    
    % we should still remove nans
%     nan_idx_final = isnan(cis_corrected);
% 
%     x(nan_idx_final) = [];
%     cis_corrected(nan_idx_final) = [];

%     cis_corrected

    figure
    plot(x, obs_interp(x))
    hold on
    plot(x, trans_interp_crxn(x))
    hold on
    plot(x, cis_corrected(x));
    hold on
    
    legend({'observed', 'trans', 'corrected'})
    title('Correcting for FPs and Transallelic Matches in FISH CDF')
    
    % save to table
    cdf_crctd_x = x;
    cdf_crctd_r = cis_corrected(x);




    cdf_corrected_table = table(cdf_crctd_x', cdf_crctd_r);
end
    
    % figure
    % plot(xq, obs_interp)
    % hold on
    % plot(cdf_obs_x, cdf_obs_f)
