 %% Finn Clark 12/6/2023
% writes a cfg file for downstream use in DNA FISH analysis

function [] = FISH_cfg_writer_for_pipe(condition_dir, file_name_root, FISH_channels, voxSize)

    basePath = condition_dir;
    baseName = file_name_root;
    projectDir = fileparts(condition_dir); %the dir containing teh dif conditions
    
    % find path to generic triclops cfg
    b = struct2table(dir(fullfile(pwd, '**/*.*')));
    t = b(contains(b.name, 'FISH_cfg_template_v2'), :);
    t_p = fullfile(t.folder{1}, t.name{1});
    temp_path = t_p;


    myLines = readlines(temp_path);

    
    for i = 1:length(myLines)
        curLine = myLines(i);
        
        % replace the placeholder for datasetpath
    
        if contains(curLine, "basePath")
    
            % replace the place holder withour files base name
            curLine = strrep(curLine, "basePath", basePath);
    
            myLines(i) = curLine;
        end

            
        if contains(curLine, "myBaseName")
    
            % replace the place holder withour files base name
            curLine = strrep(curLine, "myBaseName", baseName);
    
            myLines(i) = curLine;
        
            
        end

        if contains(curLine, "dxy_place")

            % replace the place holder withour files base name
            curLine = strrep(curLine, "dxy_place", string(voxSize(1)));
    
            myLines(i) = curLine;
        end

        if contains(curLine, "dz_place")

            % replace the place holder withour files base name
            curLine = strrep(curLine, "dz_place", string(voxSize(3)));
    
            myLines(i) = curLine;
        end

        if contains(curLine, "fishChannels_place")
        
            % replace the place holder withour files base name
            curLine = strrep(curLine, "fishChannels_place", strjoin(string(FISH_channels), ","));
            
            myLines(i) = curLine;
        end
        if contains(curLine, "refCh_place")
        
            % replace the place holder withour files base name
            curLine = strrep(curLine, "refCh_place", string(ref_channel));
            
            myLines(i) = curLine;
        end
    
    
    
    end
    
%     myLines
    
    cfg_savePath = fullfile(projectDir, baseName + ".ini");
    
    writelines(myLines, cfg_savePath);
    
    
    disp('');
    
    disp('cfg file saved to '+ string(cfg_savePath))
    
    disp('');
    
    disp('~~~~~~~');
    disp('Process Complete');
end