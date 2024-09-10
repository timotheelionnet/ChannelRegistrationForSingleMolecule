%% build the fileSet object

cfg = IniConfig();
cfg.ReadFile(config_fname);

%% Identify spots within matching distance in all pair combinations of channels


% extract the list of channels from the config object
fishChannelList = cfg.GetValues('[channelDescription]', 'fishChannelList');

outFolder = fullfile(myParams.Output.outFolder, 'res', 'matched_spots');
subOutFolder = fullfile(outFolder, 'matched_spot_loc');

% useFinnOutputFiles = cfg.GetValues('[Settings]', 'useFinnOutputFiles');

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

if ~exist(subOutFolder,'dir')
    mkdir(subOutFolder);
end
%%
nc = numel(fishChannelList); % number of channels
fovList = unique(fs.fList.FOV); % number of FOVs

% collect data on spots per unique cell and channels per cell
% cell ij has three cols: number of cells with a spot in ch i, same for j,
% number with both
sect_mat = cell(max(fishChannelList));

% init matrix to collect transallelic matching data
trans_matches = cell(1, max(fishChannelList)); % the ith element is for fish chnnel i

% loop through channels 
for i=1:nc
   % loop through all unique channel pair combinations
    for j=i+1:nc

        % will store trans matching data for channel i
        cur_trans_matches = [];
        cur_trans_matches2 = []; % use this to get the last channel
        
        % loop through fovs
        for k=1:numel(fovList)
            disp(' ');
            fileNotFound = 0;

            % load spots for channel 1

            f1 = fs.getFileName({'Channel','FOV'},...
                    {fishChannelList(i),fovList(k)},file2match,'first');
            
            [fileNotFound,loci] = loadFile(f1,'1',fileNotFound);


            
            % load spots for channel 2
            f2 = fs.getFileName({'Channel','FOV'},...
                {fishChannelList(j),fovList(k)},file2match,'first');
            
            [fileNotFound,locj] = loadFile(f2,'2',fileNotFound);


        
           % see how often we have spots from each channel in a cell
           cellID_col_indx = 7; % cell IDs are appended as the 7th column (following 6 columns of loc file)

             if ~(isempty(loci))
                 cellIDs1 = unique(loci(:,cellID_col_indx));
                 nCelli = numel(cellIDs1);
             else
                 nCelli = 0;
             end
             if ~(isempty(locj))
                 cellIDs2 = unique(locj(:,cellID_col_indx));
                 nCellj = numel(cellIDs2);
             else
                 nCellj =0;
             end
            if ~isempty(loci) && ~ isempty(locj)
                ch_intersect = intersect(cellIDs1, cellIDs2);
                % how many unique cells have both channels
                nSameCell = numel(ch_intersect);
            else
                nSameCell = 0;
            end
            sect_mat{fishChannelList(i), fishChannelList(j)} =...
                    vertcat(sect_mat{fishChannelList(i), fishChannelList(j)}, [nCelli, nCellj, nSameCell]);
            
         
            % process spots
            if ~fileNotFound
                
                if isempty(loci) || isempty (locj)
                    
                % save empty list of matched spot pairs in file 
                    t = table();
                    pairFileName = buildPairFileName(fovList,k,...
                        fishChannelList,i,j,matchDistThreshold);
                    writetable(t,fullfile(subOutFolder,pairFileName),...
                        'Delimiter','\t');
                
                    continue
                    
                end
                
                
                loci = convert_loc_pix_to_nm(loci, voxSize);
                locj = convert_loc_pix_to_nm(locj, voxSize);

                if do_cdf_crxn ==1
                    % extract transallelic matching data
                     cur_trans_matches =  [cur_trans_matches; get_trans_cdf(loci, matchDistThreshold, mutualNN_only, nDims)];
                    if j == nc
                        cur_trans_matches2 = [cur_trans_matches2; get_trans_cdf(locj, matchDistThreshold, mutualNN_only, nDims)];
                    end
                end

                
                % build distance matrix for all spots: loci vs locj
                [dr,dxyz] = compute_distMatrix_from_loc(loci,locj,nDims);
                
                % if true, we will only take spot pairs if they are mutual nearest neighbors
                if mutualNN_only ==1 

                    % get closest spot for each spot in loci
                    [min_i_vals, min_i_idx] = min(dr, [], 2);
                    
                    % get closest spot and idx for each column
                    [min_j_vals, min_j_idx] = min(dr);
                    
                    % gets the linear indices of all minima
                    mins_linearInd = sub2ind(size(dr), 1:height(dr), min_i_idx');
                    
                    % which minima are below our thresh
                    is_below_thresh = (min_i_vals < matchDistThreshold)';
                    
                    % find where min for pair ij is mutual with i as ref
                    mutual_nn_idcs  = (min_i_vals' == min(dr(:,min_i_idx))) == 1;
                    
                    keep_minima_idx = is_below_thresh & mutual_nn_idcs;
                    
                    
                    % keep only the NN ones
                    
                    matIdx = (mins_linearInd(keep_minima_idx))';
                else
    
                    % collect index pairs for all spot distances under the threshold
                    matIdx = find(dr<matchDistThreshold);
                end
                

 


                % error bypass if no spots pass 
                if isempty(matIdx) || size(matIdx,1) ==1
                   
                    t = table();
                    pairFileName = buildPairFileName(fovList,k,...
                        fishChannelList,i,j,matchDistThreshold);
                    writetable(t,fullfile(subOutFolder,pairFileName),...
                        'Delimiter','\t');
                
                    continue
                    
                end

                [r,c] = ind2sub(size(dr),matIdx);

                disp(['   C',num2str(fishChannelList(i)),' has ',num2str(size(loci,1)),...
                    ' spots; ',...
                    ' C',num2str(fishChannelList(j)),' has ',num2str(size(locj,1)),...
                    ' spots; ',...
                    'Found ',num2str(numel(r)),' spot pairs under the ',...
                    num2str(matchDistThreshold),' nm threshold.']);

                % collect coordinates of the matched spot pairs
                t = combineDistancesIntoTable(i,j,r,c,...
                    loci,locj,dr,matIdx,dxyz,fishChannelList);

                % save list of matched spot pairs in file 
                % (all spots matched regardless of uniqueness)
                pairFileName = buildPairFileName(fovList,k,...
                    fishChannelList,i,j,matchDistThreshold);
                writetable(t,fullfile(subOutFolder,pairFileName),...
                    'Delimiter','\t');
            end
        end
    end

    if do_cdf_crxn ==1
        if i<nc
            trans_matches{fishChannelList(i)} = cur_trans_matches;
        end
    end
end

if do_cdf_crxn ==1
    trans_matches{fishChannelList(nc)} = cur_trans_matches2;
end

%% generate graph of connected spots across all channels for each FOV
% channel mutliplier in order to generate a unique index for each spot
% that includes its channel of origin. e.g. channel 1 spot 8 with
% multiplier 100000 will become 100008

nc = numel(fishChannelList); % number of channels
fovList = unique(fs.fList.FOV); % list of FOVs


channelMultiplier = 100000;
fovMultiplier = 10000000;
globT = [];
for k = 1:numel(fovList) 
    for i = 1:nc  
        for j = i+1:nc   
            % load file holding the paired spots for channels i and j
            pairFileName = buildPairFileName(fovList,k,...
                    fishChannelList,i,j,matchDistThreshold);
            t = readtable(fullfile(subOutFolder,pairFileName));
            
            if isempty(t)
                
                continue
            end

            % apply channel multiplier to generate index per spot that will
            % be unique across all channels, add columns for channel and FOV or origin 
            t = reformatSpotPairTableToCombine(t,...
                fishChannelList,i,j,fovList,k,channelMultiplier,fovMultiplier);
            
            globT = [globT;t];           
        end
    end
end
% generate ranked indices to make graph generation easier and
% reformat the table of all spot pairs within all FOVs.
globT = reformatFovSpotPairTable(globT,channelMultiplier,fovMultiplier);

% generate graph object from the spot pair list
g = generateGraphFromFovTable(globT);

% graph statistics
[spotGroups,groupSizes] = conncomp(g);
[groupIDs,pairGroups,nChannelsPerGroup] = ...
    displayGroupStats(spotGroups,groupSizes,globT);

% save table of all spot pairs across the dataset, including the info about
% the ID of the connected group each pair belongs to (spotGroupID) and the
% number of spots that form the group (spotGroupSize).
globT = renamevars(globT,{'ranki','rankj'},{'spotRanki','spotRankj'});
globT = addvars(globT,pairGroups,'NewVariableNames',{'spotGroupID'});
globT = addvars(globT,groupSizes(pairGroups)','NewVariableNames',{'spotGroupSize'});
writetable(globT,fullfile(outFolder,'globalSpotPairs.txt'),'Delimiter','\t');

%% Pull out distance distribution stats
% no correction at this stage - raw data 

% all spots, no offset correction
nc = numel(fishChannelList); % number of channels
fovList = unique(fs.fList.FOV); % number of FOVs
binSize = cfg.GetValues('[Settings]', 'hist_binSize');
tmpDistArr = [];
for i = 1:nc  
        for j = i+1:nc 
            curT = globT(globT.Ci == fishChannelList(i) ...
                & globT.Cj == fishChannelList(j),:);

            % display spot pair count stats in command line
            disp(['Channel ',num2str(fishChannelList(i)),' Channel ',...
                num2str(fishChannelList(j)),'; Combined ',num2str(size(curT,1)),...
                ' spot pairs collected across ',num2str(numel(unique(curT.fovID))),...
                ' FOVs; threshold is ',num2str(matchDistThreshold),' nm']);

            % display x,y,z, offset stats
            drMed = median(curT.drij_nm);
            dxMed = median(curT.dxij_nm);
            dyMed = median(curT.dyij_nm);
            dzMed = median(curT.dzij_nm);
            tmpDistArr = [tmpDistArr;...
                fishChannelList(i),fishChannelList(j),...
                round(drMed),round(dxMed),round(dyMed),round(dzMed)];

            disp(['    Median Values: dr = ',num2str(drMed),' nm; ',...
                'dx = ',num2str(dxMed),' nm; ',...
                'dy = ',num2str(dyMed),' nm; ',...
                'dz = ',num2str(dzMed),' nm; ',...
                ]);

            figName = ['Spot distance, channels ',...
                num2str(fishChannelList(i)),'-',num2str(fishChannelList(j)),...
                'Threshold ',num2str(matchDistThreshold),' nm'];
            plotPairDistCdfs(curT,figName);
        end
end

%% add average offset to each channel to correct residual

xyz = {'x','y','z'};
distNorm = 1; 
offsetVals = zeros(numel(fishChannelList),numel(xyz));
Tcorr = globT;

if subtractAvgOffset == true
    disp(' ');
    disp('Correcting for average offset in dx dy dz');
    disp(' ');
    for i=1:numel(xyz)
        range = 1000; %start with 1 um search, then reduce range size iteratively
        while range > 0.1
            nBins = 5;
            [Tcorr, oV] = ...
                findOptimalOffset(Tcorr,xyz{i},fishChannelList,range, nBins,distNorm);
            offsetVals(:,i) = offsetVals(:,i) + oV';
            range = range/4;
        end
    end
end

if subtractAvgOffset == true
    % save table of offset-corrected coordinates and distances
    writetable(Tcorr,fullfile(outFolder,'globalSpotPairsOffsetCorrected.txt'),...
        'Delimiter','\t');
    
    % save offsets
    offsetVals = table(reshape(fishChannelList,numel(fishChannelList),1),...
        offsetVals(:,1),offsetVals(:,2),offsetVals(:,3),...
        'VariableNames',{'Channel','xOffset','yOffset','zOffset'});
    writetable(offsetVals,fullfile(outFolder,'offsetValuesPerChannel.txt'),...
        'Delimiter','\t');
end
%% Pull out distance distribution stats, after gating for spot connectivity and offset correction.
% spots in interactions, n >=4, offset correction

spotsInCell = 1;

% subtractAvgOffset = 0; % set to zero for raw coordinates; 
    % set to 1 for channel-corrected coordinates (Tcorr in previous section)
    % set to 2 for pair-corrected distances so that all median offsets
    % are strictly zero.

% plot all spots in 3D in order to spot any obvious artefacts in corners
% that could be fixed by cropping.
figure('Name','All spots xyz coordinates.');
if subtractAvgOffset == 1
    xAll =  [Tcorr.xi_nm;Tcorr.xj_nm];
    yAll =  [Tcorr.yi_nm;Tcorr.yj_nm];
    zAll =  [Tcorr.zi_nm;Tcorr.zj_nm];
else
    xAll =  [globT.xi_nm;globT.xj_nm];
    yAll =  [globT.yi_nm;globT.yj_nm];
    zAll =  [globT.zi_nm;globT.zj_nm];
end
scatter3(xAll,yAll,zAll);


%%
%cropWindowInNm = [10000, 127000; 10000, 127000];
cropWindowInNm = [0, Inf; 0, Inf]; % no cropping
%cropWindowInNm = [40000, 100000; 40000, 100000]; % super stringent window


if subtractAvgOffset == 1
    sortedT = gateSpotsOnConnectivity(...
        Tcorr,groupIDs,groupSizes,pairGroups,nChannelsPerGroup, ...
        minSpotsPerGroup,maxSpotsPerGroup,...
        minChannelsPerGroup,maxChannelsPerGroup,...
        spotsInCell,strictlyOneSpotPerChannel,cropWindowInNm);
else
    sortedT = gateSpotsOnConnectivity(...
        globT,groupIDs,groupSizes,pairGroups,nChannelsPerGroup, ...
        minSpotsPerGroup,maxSpotsPerGroup,...
        minChannelsPerGroup,maxChannelsPerGroup,...
        spotsInCell,strictlyOneSpotPerChannel,cropWindowInNm);
end

% save the table that has been gated on connectivitiy

connect_options = strcat( string(minSpotsPerGroup), '-', string(maxSpotsPerGroup), 'sptsPerGrp_', ...
                          string(minChannelsPerGroup), '-', string(maxChannelsPerGroup), 'chPerGrp' );

sortedT_fname = strcat('connectivity_gated_globalSpotPairs', connect_options, '.txt');

writetable(globT,fullfile(outFolder, sortedT_fname ),'Delimiter','\t');

disp('~~~~~~~~~~')
disp('saved connectivity gated table to')
disp(fullfile(outFolder, sortedT_fname ))
disp('~~~~~~~~~~')



nc = numel(fishChannelList); % number of channels
fovList = unique(fs.fList.FOV); % number of FOVs
binSize = cfg.GetValues('[Settings]', 'hist_binSize');
tmpDistArr = [];
for i = 1:nc  
    for j = i+1:nc 
        curT = sortedT(sortedT.Ci == fishChannelList(i) ...
            & sortedT.Cj == fishChannelList(j),:);

        % correct for median offset in x y and z and recompute dr
        if subtractAvgOffset == 2
            curT.dxij_nm = curT.dxij_nm - median(curT.dxij_nm);
            curT.dyij_nm = curT.dyij_nm - median(curT.dyij_nm);
            curT.dzij_nm = curT.dzij_nm - median(curT.dzij_nm);
            curT.drij_nm = sqrt(curT.dxij_nm.^2 + curT.dyij_nm.^2 + curT.dzij_nm.^2);
        end

        % display spot pair count stats in command line
        disp(['Channel ',num2str(fishChannelList(i)),' Channel ',...
            num2str(fishChannelList(j)),'; Combined ',num2str(size(curT,1)),...
            ' spot pairs collected across ',num2str(numel(unique(curT.fovID))),...
            ' FOVs; threshold is ',num2str(matchDistThreshold),' nm; ',...
            'Minimum spots per group: ',num2str(minSpotsPerGroup)]);

        % display x,y,z, offset stats
        drMed = median(curT.drij_nm);
        dxMed = median(curT.dxij_nm);
        dyMed = median(curT.dyij_nm);
        dzMed = median(curT.dzij_nm);
        tmpDistArr = [tmpDistArr;...
            fishChannelList(i),fishChannelList(j),...
            round(drMed),round(dxMed),round(dyMed),round(dzMed)];
    
            disp(['    Median Values: dr = ',num2str(drMed),' nm; ',...
                'dx = ',num2str(dxMed),' nm; ',...
                'dy = ',num2str(dyMed),' nm; ',...
                'dz = ',num2str(dzMed),' nm; ',...
                ]);
        if subtractAvgOffset == 1
            offsetString = ' offset corrected 1';
        elseif subtractAvgOffset == 2
            offsetString = ' offset corrected 2';
        elseif subtractAvgOffset == 0
            offsetString = ' offset corrected 0';
        end

        figName = ['Spot distance, channels ',...
            num2str(fishChannelList(i)),'-',num2str(fishChannelList(j)),...
            'Threshold ',num2str(matchDistThreshold),...
            ' nm,  minSpotsPerGroup:',num2str(minSpotsPerGroup),'; '...
            'maxSpotsPerGroup:',num2str(maxSpotsPerGroup),'; '...
            'minChannelsPerGroup: ',num2str(minChannelsPerGroup),'; ',...
            'maxChannelsPerGroup: ',num2str(minChannelsPerGroup),'; ',...
            offsetString];

        % plot and save cdfs
        cdfTable = plotPairDistCdfs(curT,figName);
        
        optionString = '';
        if strictlyOneSpotPerChannel optionString = [optionString,'_1spPerCh']; end
        if spotsInCell optionString = [optionString,'_inCellOnly']; end
        if subtractAvgOffset == 1
            optionString = [optionString,'_subtractOffset1']; 
        elseif subtractAvgOffset == 2
            optionString = [optionString,'_subtractOffset2']; 
        else
            optionString = [optionString,'_noOffsectCorxn']; 
            
        end
        
        % save tables for easy cdf plotting (saves cdf points, not raw
        % data)
        cdf_table_path = fullfile(outFolder,...
            ['cdfs',num2str(fishChannelList(i)),num2str(fishChannelList(j)),...
            '_Sp',num2str(minSpotsPerGroup),'to',num2str(maxSpotsPerGroup),...
            '_Ch',num2str(minChannelsPerGroup),'to',num2str(maxChannelsPerGroup),...
            optionString,'.txt']);
        
        writetable(cdfTable, cdf_table_path,'Delimiter','\t');


        % we want the option to do a cdf_correction
         cdf_crxn_table_path = fullfile(outFolder,...
            ['cdfs_FP-CRXN',num2str(fishChannelList(i)),num2str(fishChannelList(j)),...
            '_Sp',num2str(minSpotsPerGroup),'to',num2str(maxSpotsPerGroup),...
            '_Ch',num2str(minChannelsPerGroup),'to',num2str(maxChannelsPerGroup),...
            optionString,'.txt']);


        if do_cdf_crxn ==1
            cdf_correctedT = cdf_purifier(readtable(cdf_table_path), trans_matches{fishChannelList(i)});
            
            writetable(cdf_correctedT, cdf_crxn_table_path, 'Delimiter','\t');
        end


    end
end

%% write current analysis params to a text file

params_text = [sprintf("Spots per channel filter: min = %d, max = %d", spotThresh(1), spotThresh(2));...
                sprintf("Mutual NN only: %d", mutualNN_only); ...
                sprintf("FP Analysis: %d", do_cdf_crxn);...
                sprintf("minSpotsPerGroup = %d", minSpotsPerGroup) ;...
                sprintf("maxSpotsPerGroup = %d", maxSpotsPerGroup) ;...
                sprintf("minChannelsPerGroup = %d", minChannelsPerGroup );...
                sprintf("maxChannelsPerGroup = %d", maxChannelsPerGroup )];

time = datetime('now','TimeZone','local','Format', 'yyyyMMdd_HH_mm_ss'); % timecode to the second
params_fname = strcat(string(time), "_analysis_params.txt");

writelines(params_text, fullfile(outFolder, params_fname));

%%
% %% plot spots in 3D
% minSpotsPerGroup = 2;
% maxSpotsPerGroup = 4;
% minChannelsPerGroup = 2;
% maxChannelsPerGroup = 4;
% spotsInCell = 1;
% strictlyOneSpotPerChannel = 1;
% linearOrder = [2,4,5,3]; % linear order of FISH channels along genomic coordinates
% 
% % color of the spots, in linear order
% channelColors = [   0,  0,  204; ...
%                     204,0,  255; ...
%                     255,0,  0; ...
%                     255,205,0; ];
% channelColors = channelColors/255;
% 
% % name of the spots, in linear order
% dispName = {'Upstream Enhancer','PPARG promoter',...
%     'Downstream Enhancer','Downstream Control'};
% 
% % settings for 3D plot
% sphereRadius = 40; % radius of sphere representing each probe set in nm
% [xS,yS,zS] = sphere(20);
% cylRadius = 20; % radius of cylinder representing each link b/w probe sets in nm
% [xC,yC,zC] = cylinder(cylRadius,20);
% cylColor = 0.6+zeros(size(xC,1),size(xC,2),3); % set color to gray
% 
% % viewing angles for 3D plots
% caz = 52; 
% cel = 22;
% 
% sortedT = gateSpotsOnConnectivity(...
%     Tcorr,groupIDs,groupSizes,pairGroups,nChannelsPerGroup, ...
%     minSpotsPerGroup,maxSpotsPerGroup,...
%     minChannelsPerGroup,maxChannelsPerGroup,...
%     spotsInCell,strictlyOneSpotPerChannel,cropWindowInNm);
% 
% nSpotsPerGroup = 4;
% idxToPlot = 1:5; % indices of spot groups to plots
% sortedT = sortedT(sortedT.spotGroupSize == nSpotsPerGroup,:);
% sortedIDs = unique(sortedT.spotGroupID);
% 
% hicColorMap = zeros(256,3);
% hicColorMap(:,1) = 1;
% hicColorMap(1:end-1,2) = (0:254)/254;
% hicColorMap(1:end-1,3) = (0:254)/254;
% hicColorMap(end,:) = 0.7;
% 
% nc = numel(fishChannelList);
% hicMatrix = zeros(nc);
% for i=1:nc
%     hicMatrix(i,i) = Inf;
%     for j=i+1:nc
%         if linearOrder(i) > linearOrder(j)
%             hicMatrix(j,i) = median(sortedT.drij_nm( ...
%                 sortedT.Ci == linearOrder(j) ...
%                 & sortedT.Cj == linearOrder(i)));
%         else
%             hicMatrix(j,i) = median(sortedT.drij_nm( ...
%                 sortedT.Ci == linearOrder(i) ...
%                 & sortedT.Cj == linearOrder(j)));
%         end
%         hicMatrix(i,j) = hicMatrix(j,i);
%     end
% end
% 
% figure('Name','Median Distance Matrix'); 
% imagesc(hicMatrix);
% colormap(hicColorMap);
% clim([min(hicMatrix(~isinf(hicMatrix))),...
%     1.01*max(hicMatrix(~isinf(hicMatrix)))]);
% 
% %%
% 
% for i=1:numel(idxToPlot)
% %for i=1:1
%     gID = sortedIDs(i);
%     tmpT = sortedT(sortedT.spotGroupID == gID,:);
%     fh = figure('Name',['Spot group ',num2str(gID)]); 
%     aSpots = axes(fh,'Position',[0.1,.3,.8,.7]); hold;
%     aHiC = axes(fh,'Position',[0.05,.05,.25,.25]);
%     % get spots in each of the channels, and order them by increasing genomic
%     % coordinate
%     xyzC = getSpotCoordinatesFromGroup(tmpT,fishChannelList,linearOrder);
% 
%     refSpot = mean(xyzC(:,1:3)); % view is centered on the center of mass
%     xyzC(:,1:3) = xyzC(:,1:3) - repmat(refSpot,size(xyzC,1),1);
%     
%     % plot links between fish spots as gray cylinders
%     for j=2:size(xyzC,1)
%         dx = xyzC(j,1) - xyzC(j-1,1);
%         dy = xyzC(j,2) - xyzC(j-1,2);
%         dz = xyzC(j,3) - xyzC(j-1,3);
%         rho = sqrt( dx^2 + dy^2 + dz^2 ); %link length
%         phi = atan(dy/dx);
%         if(dx<0)
%             phi = phi+pi;
%         end
%         theta = acos(dz/rho);
%         
%         % scale the cylinder to the length of the molecular link
%         curZC0 = zC*rho;
%         curXC0 = xC;
%         curYC0 = yC;
% 
%         % rotate the cylinder by angle theta around the y axis
%         curZC1 = curZC0*cos(theta) - curXC0*sin(theta);
%         curXC1 = curZC0*sin(theta) + curXC0*cos(theta);
%         curYC1 = curYC0;
% 
%         % rotate the cylinder by phi around z axis 
%         curXC = curXC1*cos(phi) - curYC1*sin(phi);
%         curYC = curXC1*sin(phi) + curYC1*cos(phi);
%         curZC = curZC1;
% 
%         % translate the cylinder so it starts at point j-1
%         curXC = curXC + xyzC(j-1,1);
%         curYC = curYC + xyzC(j-1,2);
%         curZC = curZC + xyzC(j-1,3);
%         
%         % plot as surface
%         surf(aSpots,curXC,curYC,curZC,...
%             cylColor,'LineStyle','none',...
%             'FaceLighting','gouraud');
%     end
%     
%     % plot fish spots as colored spheres
%     for j=1:size(xyzC,1)
%         sphereColor = zeros(size(xS,1),size(xS,2),3);
%         sphereColor(:,:,1) = channelColors(j,1);
%         sphereColor(:,:,2) = channelColors(j,2);
%         sphereColor(:,:,3) = channelColors(j,3);
%         surf(aSpots,sphereRadius*xS+xyzC(j,1),sphereRadius*yS+xyzC(j,2),sphereRadius*zS+xyzC(j,3),...
%             sphereColor,'LineStyle','none','DisplayName',dispName{j},...
%             'FaceLighting','gouraud');
%     end
% 
%     % set plot legends, axes etc.
%     light(aSpots);
%     grid(aSpots,'on');
%     xlabel(aSpots,'x in nm');
%     ylabel(aSpots,'y in nm');
%     zlabel(aSpots,'z in nm');
%     axis(aSpots,'equal');
%     xlim(aSpots,[-600,600]);
%     ylim(aSpots,[-600,600]);
%     zlim(aSpots,[-200,200]);
%     view(aSpots,[caz,cel]);
% 
%     % plot HiC-like distance matrix
%     curHicMatrix = hicMatrix;
%     for j=1:nc
%         for k=j+1:nc
%             curHicMatrix(j,k) = sqrt(sum( (xyzC(xyzC(:,4) == linearOrder(j) ,1:3) ...
%                     - xyzC(xyzC(:,4) == linearOrder(k) ,1:3)).^2 ));
%         end
%     end
%     
%     cMin = 200;
%     cMax = 400;
%     imagesc(aHiC,cMax - (cMax-cMin)/255 + zeros(4,4)); 
%     hold;
%     h = imagesc(aHiC,curHicMatrix);
%     alphaMatrix = ones(4);
%     alphaMatrix(curHicMatrix > cMax) = 0;
%     alphaMatrix(logical(eye(4))) = 1;
%     set(h,'AlphaData',alphaMatrix);
%     colormap(hicColorMap);
% %     clim([min(curHicMatrix(~isinf(curHicMatrix))),...
% %         1.01*max(curHicMatrix(~isinf(curHicMatrix)))]);
% %     clim([min(hicMatrix(~isinf(hicMatrix))),...
% %         1.01*max(hicMatrix(~isinf(hicMatrix)))]);
%      clim([cMin,cMax]);
%     colorbar;
% end

%clear curXC0 curYC0 curZC0 curXC1 curYC1 curZC1 caz cel xS yS zS xC yC zC;

%% ~~~~~~~~~~~~~~~end of main script body~~~~~~~~~~~~~~~~~
% below are subfunctions

%% 
function xyzC = getSpotCoordinatesFromGroup(tmpT,fishChannelList,linearOrder)
    xyz = zeros(numel(fishChannelList),3);
    for i=1:numel(fishChannelList)
        idxi = find(ismember(tmpT.Ci,fishChannelList(i)),1);
        if ~isempty(idxi)   
            xyz(i,:) = [tmpT.xi_nm(idxi),tmpT.yi_nm(idxi),tmpT.zi_nm(idxi)];
        else
            idxj = find(ismember(tmpT.Cj,fishChannelList(i)),1) ;
            if ~isempty(idxj)
                xyz(i,:) = [tmpT.xj_nm(idxj),tmpT.yj_nm(idxj),tmpT.zj_nm(idxj)];
            end
        end
    end
    xyz = [xyz,fishChannelList'];
    xyzC = zeros(size(xyz));
    for i=1:numel(linearOrder)
        xyzC(i,:) = xyz(xyz(:,4)==linearOrder(i),:);
    end
end


%% find offset coordinate that minimizes the pair distances across channels
% parameter sweep to find the offsets that when added to the different channels 
% minimize the total error in the pair distances.
% T is input table of coordinates and distances
% xyz = 'x' or 'y' or 'z' sets which coordinate to offset
% range is the range in nm of the parameter sweep (will sweep values from -range to range)
% nBins is the half number of bins
% distNorm is the distance norm used for minimization 
    % set to 1 for sum of absolute median distances across all pairs of
    % channels
    % set to 2 for sum of median distance squares across all pairs of
    % channels
function [Tcorr, offsetVals] = findOptimalOffset(...
                T,xyz,fishChannelList,range, nBins,distNorm) 
    nBins = 2*nBins;
    disp([xyz,' coordinate; testing offsets from ',num2str(-range),...
        ' nm to ',num2str(range),' nm with ',...
        num2str((2*range)/nBins),' nm steps...']);
    xOffset = -range:(2*range)/nBins:range;
    
    if numel(fishChannelList) == 4
        [dxMin0,dyMin0,dzMin0] =  compute_dxdydz(T,fishChannelList,distNorm);
        switch xyz
            case 'x'
                dxyzMin0 = dxMin0;
            case 'y'
                dxyzMin0 = dyMin0;
            case 'z'
                dxyzMin0 = dzMin0;
        end
        iMin3 = nBins/2+1;
        iMin4 = nBins/2+1;
        iMin5 = nBins/2+1;
        dxyzMin = dxyzMin0;
        for i3= 1:numel(xOffset) 
            for i4=1:numel(xOffset)
                for i5=1:numel(xOffset)    
                    offsetVals = [0,xOffset(i3),xOffset(i4),xOffset(i5)];
                    cT = addOffsetToXYZCoordinates(...
                        T,offsetVals,xyz,fishChannelList);
                    [dx,dy,dz] = compute_dxdydz(cT,fishChannelList,distNorm);
                    switch xyz
                        case 'x'
                            dxyz = dx;
                        case 'y'
                            dxyz = dy;
                        case 'z'
                            dxyz = dz;
                    end
                    if(dxyz < dxyzMin)
                        dxyzMin = dxyz;
                        iMin3 = i3;
                        iMin4 = i4;
                        iMin5 = i5;
                    end
                end
            end
        end
        offsetVals = [0,xOffset(iMin3),xOffset(iMin4),xOffset(iMin5)];
        Tcorr = addOffsetToXYZCoordinates(...
                        T,offsetVals,xyz,fishChannelList);
        disp(['initial avg d',xyz,' = ',num2str(dxyzMin0),... 
            '; final d',xyz,' = ',num2str(dxyzMin),...
            '; offset',num2str(fishChannelList(2)),' = ',...
            num2str(offsetVals(2)),' nm; ',...
            '; offset',num2str(fishChannelList(3)),' = ',...
            num2str(offsetVals(3)),' nm; ',...
            '; offset',num2str(fishChannelList(3)),' = ',...
            num2str(offsetVals(4)),' nm; ']);      
    end
end


%% add offset in different channels
% from an input table Tin, 
% this functions adjusts the coordinate xyz (enter either 'x', 'y' or 'z')
% by an offset given by offsetVals array, which should have the same number
% of elements as fishChannelList.
% e.g. if xyz = 'z',  offsetVals = [0,0,100,0] and fishChannelList = [2,3,4,5]
% the function adds an offset of 100 nm to the z coordinate of channel 4,
% it updates the pair distances in the table in the process.
function Tout = addOffsetToXYZCoordinates(Tin,offsetVals,xyz,fishChannelList)
    % find variable name for the requested coordinate
    curXYZi = [xyz,'i_nm'];
    curXYZj = [xyz,'j_nm'];

    % adjust coordinate value by offset
    Tout = Tin;
    nc = numel(fishChannelList);
    for i =1:nc
        Tout.(curXYZi)(Tin.Ci==fishChannelList(i)) = ...
            Tin.(curXYZi)(Tin.Ci==fishChannelList(i)) + offsetVals(i);
        Tout.(curXYZj)(Tin.Cj==fishChannelList(i)) = ...
            Tin.(curXYZj)(Tin.Cj==fishChannelList(i)) + offsetVals(i);
    end

    % once coordinates have been adjusted for offset, make sure the
    % distances are up to date.
    deltaXYZ = ['d',xyz,'ij_nm'];
    Tout.(deltaXYZ) = Tout.(curXYZi) - Tout.(curXYZj);
    Tout.drij_nm = ...
        sqrt(Tout.dxij_nm.^2 + Tout.dyij_nm.^2 + Tout.dzij_nm.^2);

end

%% compute a metric of the total distances across channel pairs
% T is the table of coordinates and pair distances
% distNorm is the distance norm used for minimization 
    % set to 1 for sum of absolute median distances across all pairs of
    % channels
    % set to 2 for sum of median distance squares across all pairs of
    % channels 
function [dx,dy,dz] = compute_dxdydz(T,fishChannelList,distNorm)
    nc = numel(fishChannelList);
    dx = 0; dy = 0; dz = 0;
    for i=1:nc
        for j=i+1:nc
            idx = (T.Ci == fishChannelList(i)) & (T.Cj == fishChannelList(j));
            if distNorm == 2
                dx = dx + median(T.xj_nm(idx)-T.xi_nm(idx))^2;
                dy = dy + median(T.yj_nm(idx)-T.yi_nm(idx))^2;
                dz = dz + median(T.zj_nm(idx)-T.zi_nm(idx))^2;
            elseif distNorm == 1
                dx = dx + abs(median(T.xj_nm(idx)-T.xi_nm(idx)));
                dy = dy + abs(median(T.yj_nm(idx)-T.yi_nm(idx)));
                dz = dz + abs(median(T.zj_nm(idx)-T.zi_nm(idx)));
            end
        end
    end
    if distNorm == 2
        dx = sqrt(dx);
        dy = sqrt(dy);
        dz = sqrt(dz);
    end
end

%% combines all metrics of spot pairs between two channels into a single table
% variables ordered as follows:
% col 1: idx<i> (e.g. idx1)
% col 2: idx<j> 
% col 3: dr<i><j>_nm
% col 4: dx<i><j>_nm
% col 5: dy<i><j>_nm
% col 6: dz<i><j>_nm
% col 7: x<i>_nm
% col 8: y<i>_nm
% col 9: z<i>_nm
% col 10: x<j>_nm
% col 11: y<j>_nm
% col 12: z<j>_nm
% col 13: cellID
function t = combineDistancesIntoTable(i,j,r,c,loci,locj,dr,matIdx,dxyz,fishChannelList)
    xi = loci(r,1);
    yi = loci(r,2);
    zi = loci(r,3);

    xj = locj(c,1);
    yj = locj(c,2);
    zj = locj(c,3);

    cellID = loci(r,7);
    
    % collect distance between spots in the pair
    drij = dr(matIdx);
    
    drij = reshape(drij,size(drij,1),1);
    
    % collect deltaX,Y,Z between spots in the pair
    dxij = dxyz{1}(matIdx);
    dxij = reshape(dxij,size(dxij,1),1);

    dyij = dxyz{2}(matIdx);
    dyij = reshape(dyij,size(dyij,1),1);
    
    dzij = dxyz{3}(matIdx);
    dzij = reshape(dzij,size(dzij,1),1);
    
    % compile all info into a table
    t = table(r,c,drij,dxij,dyij,dzij,xi,yi,zi,xj,yj,zj,cellID);
    t = renamevars(t,{'r','c'},...
        {['idx',num2str(fishChannelList(i))],...
         ['idx',num2str(fishChannelList(j))]});
    t = renamevars(t,{'xi','yi','zi'},...
        {['x',num2str(fishChannelList(i)),'_nm'],...
         ['y',num2str(fishChannelList(i)),'_nm'],...
         ['z',num2str(fishChannelList(i)),'_nm']});
    t = renamevars(t,{'xj','yj','zj'},...
        {['x',num2str(fishChannelList(j)),'_nm'],...
         ['y',num2str(fishChannelList(j)),'_nm'],...
         ['z',num2str(fishChannelList(j)),'_nm']});
    t = renamevars(t,{'drij','dxij','dyij','dzij'},...
        {['dr',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm'],...
         ['dx',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm'],...
         ['dy',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm'],...
         ['dz',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm']});
end

%% plot a histogram of how many pairs form connected groups of different sizes
% and how many channels are present in each connected groups
function [groupIDs,pairGroups,nChannelsPerGroup] = displayGroupStats(spotGroups,groupSizes,T)
    
    % count all spots
    s = unique(union(T.ranki,T.rankj));
    
    disp(['There are ',num2str(numel(s)),' unique spots in total, ',...
        num2str(size(T,1)),' spot pairs, and ',num2str(numel(groupSizes)),' groups.']);
    
    % internal check that the sum of all spot pairs per groups adds up to all
    % the pairs
    nS = sum(groupSizes);
    disp(['    Internal check: number of spot pairs recovered from group sizes: ',...
        num2str(nS)])
    
    %% plot histogram of group sizes
    fh = figure('Name',...
        ['Histogram of group sizes (',num2str(numel(groupSizes)),' total groups)']);
    [n,x] = hist(groupSizes,0:max(groupSizes)+1);
    bar(x,n);
    xlabel('Group size');
    ylabel('Number of groups');
    
    %% plot histogram of number of channels per group
    
    % find the group that each spot pair belongs to in the table T
    pairGroups = spotGroups(T.ranki)';
    
    % internal control: check that spot i index and spot j index belong to the
    % same group
    pg2 = spotGroups(T.rankj)';
    if sum((pairGroups-pg2).^2)>0
        disp('Failed check that spot i index and spot j index belong to the same group');
    end
    
    % compute histogram 
    groupIDs = unique(spotGroups);
    nChannelsPerGroup = zeros(size(groupIDs));
    for i=1:numel(groupIDs)
        nChannelsPerGroup(i) = numel(union( ...
            unique(T.Ci(pairGroups==groupIDs(i))),...
            unique(T.Cj(pairGroups==groupIDs(i))) ));
    end
    
    fh = figure('Name',['Histogram of number of channels per group (',num2str(numel(groupSizes)),' total groups)']);
    [n,x] = hist(nChannelsPerGroup,0:max(nChannelsPerGroup)+1);
    bar(x,n);
    xlabel('Number of channels per group');
    ylabel('Number of groups');

    nCC = sum(nChannelsPerGroup == groupSizes);
    disp([num2str(nCC),' Groups have striclty one spot per channel (',num2str(100*nCC/numel(groupSizes)),'%).']);

    nCC = sum(nChannelsPerGroup == groupSizes & groupSizes == 4);
    disp([num2str(nCC),' Groups have striclty one spot in all of the channels (',num2str(100*nCC/numel(groupSizes)),'%).']);

end


%% gate interactions based on group size or composition

function sortedT = gateSpotsOnConnectivity(...
    T,groupIDs,groupSizes,pairGroups,nChannelsPerGroup, ...
    minSpotsPerGroup,maxSpotsPerGroup,...
    minChannelsPerGroup,maxChannelsPerGroup,...
    spotsInCell,strictlyOneSpotPerChannel,cropWindowInNm)
    
    selectedGroupIDs = (groupSizes >= minSpotsPerGroup) ...
        & (groupSizes <= maxSpotsPerGroup);
    size(selectedGroupIDs)
    sum(selectedGroupIDs)

    selectedGroupIDs = selectedGroupIDs ...
        & (nChannelsPerGroup >= minChannelsPerGroup) ...
        & (nChannelsPerGroup <= maxChannelsPerGroup);
    sum(selectedGroupIDs)

    if strictlyOneSpotPerChannel
        selectedGroupIDs = selectedGroupIDs ...
            & ( nChannelsPerGroup == groupSizes);
    end
    sum(selectedGroupIDs)

    selectedGroupIDs = groupIDs(selectedGroupIDs);
    
    min(selectedGroupIDs)
    max(selectedGroupIDs)
    numel(selectedGroupIDs)
    
    idx = ismember(pairGroups,selectedGroupIDs);
    if spotsInCell
        spotsInCellString = 'spots in cells only';
        idx = idx & (T.cellID ~=0);
    else
        spotsInCellString = '';
    end

    nT = size(T,1);
    isWithinWindow = (T.xi_nm >= repmat(cropWindowInNm(1,1),nT,1)) ...
        & (T.xi_nm <= repmat(cropWindowInNm(1,2),nT,1)) ...
        & (T.yi_nm >= repmat(cropWindowInNm(2,1),nT,1)) ...
        & (T.yi_nm <= repmat(cropWindowInNm(2,2),nT,1));
    
    idx = idx & isWithinWindow;
    
    sortedT = T( idx,:);
    
    n = size(sortedT,1);
    disp(['Selected ',num2str(n),'/',num2str(size(T,1)),...
        ' spots based on minSpotsPerGroup = ',num2str(minSpotsPerGroup),...
        '; maxSpotsPerGroup = ',num2str(maxSpotsPerGroup),...
        '; minChannelsPerGroup = ',num2str(minChannelsPerGroup),'; ',...
        'maxChannelsPerGroup = ',num2str(maxChannelsPerGroup),'; ', ...
        spotsInCellString ,...
        'strictlyOneSpotPerChannel = ', num2str(strictlyOneSpotPerChannel),...
        'cropXmin_nm = ',num2str(cropWindowInNm(1,1)),...
        'cropXmax_nm = ',num2str(cropWindowInNm(2,1)),...
        'cropYmin_nm = ',num2str(cropWindowInNm(1,2)),...
        'cropYmax_nm = ',num2str(cropWindowInNm(2,2)),...
        ]);
end

%% plot dx dy dz cumulated distribution functions;
% also outputs a table of all cdf values for easy saving.
function cdfTable = ...
    plotPairDistCdfs(curT,figName)

    % create figure and plots
    fh = figure('Name',figName);
    
    aR = axes(fh,'Position',[0.1,0.1,0.35,0.8]); % dr between channels
    aX = axes(fh,'Position',[0.55,0.7,0.35,0.2]); % dx
    aY = axes(fh,'Position',[0.55,0.4,0.35,0.2]); % dy 
    aZ = axes(fh,'Position',[0.55,0.1,0.35,0.2]); % dz
    
    % add random epsilon to avoid distances happening to be identical
    eps = 1e-9*median(curT.drij_nm);
    curT.drij_nm = curT.drij_nm + eps*rand(size(curT,1),1);
    curT.dxij_nm = curT.dxij_nm + eps*rand(size(curT,1),1);
    curT.dyij_nm = curT.dyij_nm + eps*rand(size(curT,1),1);
    curT.dzij_nm = curT.dzij_nm + eps*rand(size(curT,1),1);

    % compute cdfs
    [cdfR_f,cdfR_x] = ecdf(curT.drij_nm);
    [cdfX_f,cdfX_x] = ecdf(curT.dxij_nm);
    [cdfY_f,cdfY_x] = ecdf(curT.dyij_nm);
    [cdfZ_f,cdfZ_x] = ecdf(curT.dzij_nm);

    
    cdfTable = table(cdfR_f,cdfR_x,cdfX_f,cdfX_x,cdfY_f,cdfY_x,cdfZ_f,cdfZ_x);

    % plot cdfs
    plot(aR,cdfR_x,cdfR_f);
    plot(aX,cdfX_x,cdfX_f);
    plot(aY,cdfY_x,cdfY_f);
    plot(aZ,cdfZ_x,cdfZ_f);
    
    % set legends and title of each graph
    aR.XGrid = 'on';
    aR.YGrid = 'on';
    aR.XLabel.String = 'Distance in nm';
    aR.YLabel.String = 'CDF';
    aR.Title.String ='dist in nm';

    aX.XGrid = 'on';
    aX.YGrid = 'on';
    aX.XLabel.String = 'Distance in nm';
    aX.YLabel.String = 'CDF';
    aX.Title.String ='X offset in nm';
    
    aY.XGrid = 'on';
    aY.YGrid = 'on';
    aY.XLabel.String = 'Distance in nm';
    aY.YLabel.String = 'CDF';
    aY.Title.String ='Y offset in nm';
   
    aZ.XGrid = 'on';
    aZ.YGrid = 'on';
    aZ.XLabel.String = 'Distance in nm';
    aZ.YLabel.String = 'CDF';
    aZ.Title.String ='Z offset in nm';
    
end



%% convert the table of spot pairs into a graph object in order to extract connected components more easily
function g = generateGraphFromFovTable(fovT)
    EndNodes = [fovT.('ranki'),fovT.('rankj')];
    EdgeTable = table(EndNodes,...
        fovT.drij_nm,fovT.dxij_nm,fovT.dyij_nm,fovT.dzij_nm,...
        fovT.Ci,fovT.Cj,fovT.fovID,...
        'VariableNames',...
        {'EndNodes','drij_nm','dxij_nm','dyij_nm','dzij_nm','Ci','Cj','fovID'});
    
    NodeTable = table(fovT.('ranki'),fovT.('idxi'),...
        fovT.('xi_nm'),fovT.('yi_nm'),fovT.('zi_nm'),...
        fovT.('Ci'),fovT.('fovID'),'VariableNames',...
        {'rank','idx','x_nm','y_nm','z_nm','C','fovID'});

    NodeTable = [NodeTable; ...
        table(fovT.('rankj'),fovT.('idxj'),...
        fovT.('xj_nm'),fovT.('yj_nm'),fovT.('zj_nm'),...
        fovT.('Cj'),fovT.('fovID'),'VariableNames',...
        {'rank','idx','x_nm','y_nm','z_nm','C','fovID'})];

    NodeTable = unique(NodeTable,'rows');

    g = graph(EdgeTable,NodeTable);
end

%% reformat the table holding the spot pair data compiled over all channel pairs of a FOV
% this changes the spot IDs which encode channel of origin into a rank
    % index ranki and rankj (simplifies the graph generation later on)
% reverts the spot index idxi idxj to their original values without the
    % channel multiplier (channel of origin should still be encoded in the variables
    % Ci and Cj)
function fovT = reformatFovSpotPairTable(fovT,channelMultiplier,fovMultiplier)
    m = min(channelMultiplier,fovMultiplier);

    fovT = renamevars(fovT,'idxi','multIdxi');
    fovT = renamevars(fovT,'idxj','multIdxj');
    idxi = mod(fovT.multIdxi,m);
    idxj = mod(fovT.multIdxj,m);
    fovT = addvars(fovT,idxi,'Before','multIdxi');
    fovT = addvars(fovT,idxj,'Before','multIdxi');
    
    % find indices of unique spots in the table (indices including the multiplier for the channel specificiation)
    allIdx = unique([fovT.('multIdxi');fovT.('multIdxj')]);

    % map spot indices to their rank order
    map = containers.Map(num2cell(allIdx),num2cell(1:numel(allIdx)));

    % replace spot indices with their mapped rank order in fovT
    ranki = values(map,num2cell(fovT.('multIdxi')));
    ranki = cell2mat(ranki);
    rankj = values(map,num2cell(fovT.('multIdxj')));
    rankj = cell2mat(rankj);
    fovT = addvars(fovT,ranki,'Before','idxi');
    fovT = addvars(fovT,rankj,'Before','idxi');
    
    fovT = removevars(fovT,'multIdxi');
    fovT = removevars(fovT,'multIdxj');
end

%% adds a channel multiplier to each spot index to encode its channel of origin and fov of origin
% this avoids spot indices ambiguities between channels.
% also removes the channel specifications in the variable names, 
% e.g. idx3 and idx4 become idxi and idxj for a channel 3/4 pair table;
% dr34_nm becomes drij_nm, etc.
% adds a channel of origin variable Ci and Cj, and a fov or origin variable
% fovID.
function t = reformatSpotPairTableToCombine(t,...
                fishChannelList,i,j,fovList,k,channelMultiplier,fovMultiplier)
    
    % reformat the spot indices so that there are no repeat indices across
    % channels e.g. channel 3 spot 8 with multiplier 100000 will become 300008
    t = applyChannelAndFovMultiplierToSpotIndices(t,...
                fishChannelList,i,j,channelMultiplier,fovList,k,fovMultiplier);

    % remove the specific channel index values from the variable names
    % e.g. idx3 and idx4 become idxi and idxj for a channel 3/4 pair table;
    % dr34_nm becomes drij_nm, etc.
    t = removeIndicesFromVariableNames(t,...
                fishChannelList,i,j);
    
    Ci = repmat(fishChannelList(i),size(t,1),1);
    Cj = repmat(fishChannelList(j),size(t,1),1);
    fovID = repmat(fovList(k),size(t,1),1);
    t = addvars(t,Ci);
    t = addvars(t,Cj);
    t = addvars(t,fovID);
    
end

%% change channel specific variable names into generic variable names
% e.g. the variable names 'idx1' and 'idx2' in table t become 'idxi' and 'idxj'
% dr12_nm becomes drij_nm, etc.
function t = removeIndicesFromVariableNames(t, fishChannelList,i,j)

    % check that variables in t have the expected names - this is needed
    % because the later part of the function expects a certain variable in a
    % certain column without checking.
    verbose= 1;
    areVariablesCorrect = checkVariableNamesInSpotPairTable(...
        t,fishChannelList,i,j,verbose);
    if ~areVariablesCorrect
        disp(['Wrong table format; ',...
            'cannot reformat table with channel-independent variable names']);
        return
    end
    
    % idx3 and idx4 become idxi and idxj for a channel 3/4 pair table;
    % dr34_nm becomes drij_nm, etc.
    t = renamevars(t,{['idx',num2str(fishChannelList(i))]},{'idxi'}); 
    t = renamevars(t,{['idx',num2str(fishChannelList(j))]},{'idxj'}); 
    t = renamevars(t,{['dr',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm']},...
        {'drij_nm'}); 
    t = renamevars(t,{['dx',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm']},...
        {'dxij_nm'}); 
    t = renamevars(t,{['dy',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm']},...
        {'dyij_nm'}); 
    t = renamevars(t,{['dz',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm']},...
        {'dzij_nm'}); 
    t = renamevars(t,{['x',num2str(fishChannelList(i)),'_nm']},{'xi_nm'}); 
    t = renamevars(t,{['y',num2str(fishChannelList(i)),'_nm']},{'yi_nm'}); 
    t = renamevars(t,{['z',num2str(fishChannelList(i)),'_nm']},{'zi_nm'}); 
    t = renamevars(t,{['x',num2str(fishChannelList(j)),'_nm']},{'xj_nm'}); 
    t = renamevars(t,{['y',num2str(fishChannelList(j)),'_nm']},{'yj_nm'}); 
    t = renamevars(t,{['z',num2str(fishChannelList(j)),'_nm']},{'zj_nm'}); 
end

%% checks that the table of Spot Pair info for channels i and j has the right format
% expects the following variable names in that order:
% col 1: idx<i> (e.g. idx1)
% col 2: idx<j> 
% col 3: dr<i><j>_nm
% col 4: dx<i><j>_nm
% col 5: dy<i><j>_nm
% col 6: dz<i><j>_nm
% col 7: x<i>_nm
% col 8: y<i>_nm
% col 9: z<i>_nm
% col 10: x<j>_nm
% col 11: y<j>_nm
% col 12: z<j>_nm
% col 13: cellID
function areVariablesCorrect = checkVariableNamesInSpotPairTable(t,fishChannelList,i,j,verbose)
    
    % this is the expected list of variables in the spot pair table
    expectedVars = {['idx',num2str(fishChannelList(i))], ...
                    ['idx',num2str(fishChannelList(j))], ...
                    ['dr',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm'],...
                    ['dx',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm'],...
                    ['dy',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm'],...
                    ['dz',num2str(fishChannelList(i)),num2str(fishChannelList(j)),'_nm'],...
                    ['x',num2str(fishChannelList(i)),'_nm'],...
                    ['y',num2str(fishChannelList(i)),'_nm'],...
                    ['z',num2str(fishChannelList(i)),'_nm'],...
                    ['x',num2str(fishChannelList(j)),'_nm'],...
                    ['y',num2str(fishChannelList(j)),'_nm'],...
                    ['z',num2str(fishChannelList(j)),'_nm'],...
                    'cellID'};

    
    % check that the expected variables are present in the right order in table t
    areVariablesCorrect = 1;
    if numel(t.Properties.VariableNames) ~= numel(expectedVars)
        areVariablesCorrect = 0;
    else
        for i=1:numel(expectedVars)
            if ~strcmp(expectedVars{i},t.Properties.VariableNames{i})
                areVariablesCorrect = 0;
            end
        end
    end

    if (~areVariablesCorrect && verbose)
        disp(['Spot Pair Table for channels ',...
            num2str(fishChannelList(i)),' and ',...
            num2str(fishChannelList(j)),...
            ' does not have the expected variables names.']);
        disp('Expected Variables:');
        disp(expectedVars);
        disp('Actual Variables:');
        disp(t.Properties.VariableNames);
    end
end

%% from a table of spot pairs, changes the indices in the first two columns 
% by adding an offset equal to the channel multiplier times the channel ID
% and the FOV multiplier times the FOV ID
% this should generate a unique index for each spot across all channels.
% e.g. channel 3 spot 8 fov 11 with 
% channel multiplier 100,000 and fov multiplier 10,000,000
% will become 110,300,008

function t = applyChannelAndFovMultiplierToSpotIndices(t,...
    fishChannelList,i,j,channelMultiplier,fovList,k,fovMultiplier)

    % check that variables in t have the expected names - this is needed
    % because the later part of the function expects a certain variable in a
    % certain column without checking.
    verbose= 1;
    areVariablesCorrect = checkVariableNamesInSpotPairTable(...
        t,fishChannelList,i,j,verbose);
    if ~areVariablesCorrect
        disp(['Wrong table format; ',...
            'cannot reformat table with channel-independent variable names']);
        return
    end

    t.(1) = t.(1) + fishChannelList(i) * channelMultiplier + fovList(k) * fovMultiplier;
    t.(2) = t.(2) + fishChannelList(j) * channelMultiplier + fovList(k) * fovMultiplier;
end

%% build pair filename
% standardized file name for spot pairs with the format
% fov<k>_C<i><j>_spotPairsMatched<matchDistThreshold>nmThresh.txt
% e.g. 
% fov1_C23_spotPairsMatched750nmThresh.txt
function pairFileName = buildPairFileName(fovList,k,fishChannelList,i,j,...
    matchDistThreshold)

pairFileName = ['fov',num2str(fovList(k)),'_'...
                    'C',num2str(fishChannelList(i)),...
                    num2str(fishChannelList(j)),'_spotPairsMatched',...
                    num2str(matchDistThreshold),'nmThresh.txt'];

end

%% load loc file
% includes safety nets if file doesnt exist
function [fileNotFound,loc] = loadFile(f,fNum,fileNotFound)
% disp(['Opening file ',fNum,': ',f{1},'...']);
    if ~isempty(f)
        if(exist(f{1},'file'))
            loc = load(f{1});
                
        else
            fileNotFound = 1;
            disp('File missing!');
            
        end
    else
        fileNotFound = 1;
        disp('File missing!');
        
    end
end

%% convert voxel units in loc file to nm
% self explanatory
function loc = convert3DLocToNm(loc,voxSize_dxy,voxSize_dz)
    loc(:,1:2) = loc(:,1:2)*voxSize_dxy;
    loc(:,3) = loc(:,3)*voxSize_dz;

end 