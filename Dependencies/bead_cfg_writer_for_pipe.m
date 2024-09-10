%% Finn Clark 12/6/2023
% writes a cfg file for downstream use in DNA FISH analysis

function [cfg_path] = bead_cfg_writer_for_pipe(code_dir, split_ch_dir, file_name_root, bead_channels, ref_channel)
    
    [d,~] = fileparts(split_ch_dir);
    basePath = split_ch_dir;
    baseName = file_name_root;
    projectDir = d; %the dir containing your project
    
    % find path to generic triclops cfg
    b = struct2table(dir(fullfile(code_dir, '**/*.*')));
    t = b(contains(b.name, 'beads_cfg_template_v2'), :);
    t_p = fullfile(t.folder{1}, t.name{1});
    temp_path = t_p;
    
    
    res_path = fullfile(d, 'res');
    
    % temp_path
    
    
    myLines = readlines(temp_path);
    
    for i = 1:length(myLines)
        curLine = myLines(i);
        
        % replace the placeholder for datasetpath
    
        if contains(curLine, "beadImg = basePath")
    
            % replace the place holder withour files base name
            curLine = strrep(curLine, "basePath", fullfile(basePath, strcat('C{Channel}-',baseName, '_BEADS-{FOV}.tif')) );
    
            myLines(i) = curLine;
        end
        
        if contains(curLine, "loc = basePath")
    
            % replace the place holder withour files base name
            curLine = strrep(curLine, "basePath", fullfile(basePath, strcat('C{Channel}-', baseName, '_BEADS-{FOV}.loc4')));
    
            myLines(i) = curLine;
        end

        
        if contains(curLine, "outFolder= basePath")
    
            % replace the place holder withour files base name
            curLine = strrep(curLine, "basePath",res_path);
    
            myLines(i) = curLine;
        end
        

        

        if contains(curLine, "fishChannels_place")
        
            % replace the place holder withour files base name
            curLine = strrep(curLine, "fishChannels_place", strjoin(string(bead_channels), ","));
            
            myLines(i) = curLine;
        end
        if contains(curLine, "refCh_place")
        
            % replace the place holder withour files base name
            curLine = strrep(curLine, "refCh_place", string(ref_channel));
            
            myLines(i) = curLine;
        end
    
    
    
    end
    
    % myLines
    
    cfg_savePath = fullfile(projectDir, baseName + ".ini");
    
    writelines(myLines, cfg_savePath); 

    cfg_path = cfg_savePath; 
    
    
    
    
    disp('');
    
    disp('cfg file saved to '+ string(cfg_savePath))
    
    disp('');
    
    disp('~~~~~~~');
    % disp('Process Complete');
end