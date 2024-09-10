%% FISH Pipeline GUI
function refinedGUI
    % Create the main figure window
    fig = uifigure('Position', [100 100 400 600], 'Name', 'Refined GUI');
    
    % Header for Section 1: Config
    header1 = uilabel(fig, 'Position', [20 570 360 20], 'Text', 'Config', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Section 1: Long text box with updated label and two small text boxes
    lbl1 = uilabel(fig, 'Position', [20 540 250 30], 'Text', 'Paste the path to your project directory:');
    longTextBox = uieditfield(fig, 'text', 'Position', [20 510 360 30], 'Value', 'path/to/your/project');
    
    lbl2a = uilabel(fig, 'Position', [20 480 150 30], 'Text', 'Voxel Size XY (nm):');
    smallTextBox1 = uieditfield(fig, 'numeric', 'Position', [20 450 150 30], 'Value', 98);
    
    lbl2b = uilabel(fig, 'Position', [230 480 150 30], 'Text', 'Voxel Size Z (nm):');
    smallTextBox2 = uieditfield(fig, 'numeric', 'Position', [230 450 150 30], 'Value', 250);
    
    % Divider 1: Horizontal line
    line1 = uipanel(fig, 'Position', [20 430 360 1], 'BackgroundColor', 'black');
    
    % Header for Section 2: Bead Correction
    header2 = uilabel(fig, 'Position', [20 420 360 20], 'Text', 'Bead Correction', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Section 2: Three checkboxes with updated labels and behavior
    chkBox1 = uicheckbox(fig, 'Position', [20 380 250 30], 'Text', 'Use Bead Buddy?', 'Value', true, 'ValueChangedFcn', @(cb,event) chkBox1Changed(cb));
    chkBox2 = uicheckbox(fig, 'Position', [20 340 250 30], 'Text', 'Do pseudo flatfield crxn?', 'Value', true, 'Enable', 'on');
    chkBox3 = uicheckbox(fig, 'Position', [20 300 250 30], 'Text', 'Perform localization?', 'Value', true, 'Enable', 'on');
    
    % Divider 2: Horizontal line
    line2 = uipanel(fig, 'Position', [20 290 360 1], 'BackgroundColor', 'black');
    
    % Header for Section 3: Image Processing
    header3 = uilabel(fig, 'Position', [20 280 360 20], 'Text', 'Image Processing', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Section 3: Two checkable boxes with updated labels and pre-filled text boxes
    chkBox4 = uicheckbox(fig, 'Position', [20 240 180 30], 'Text', 'Make cell masks?', 'Value', true, 'ValueChangedFcn', @(cb,event) chkBox4Changed(cb));
    chkBox5 = uicheckbox(fig, 'Position', [20 200 180 30], 'Text', 'Do pseudo flatfield crxn?', 'Value', true);
    
    lbl4 = uilabel(fig, 'Position', [20 170 120 30], 'Text', 'Mask Channel:');
    maskChannel = uieditfield(fig, 'numeric', 'Position', [20 140 150 30], 'Value', 1);
    
    lbl5 = uilabel(fig, 'Position', [230 170 180 30], 'Text', 'Average cell diam (pixels):');
    avgCellDiam = uieditfield(fig, 'numeric', 'Position', [230 140 150 30], 'Value', 150);
    
    % Save and Run button
    btnSaveRun = uibutton(fig, 'push', 'Position', [120 20 160 30], 'Text', 'Save and Run', 'ButtonPushedFcn', @(btn,event) saveAndRun());
    
    % Pause script and wait for user interaction
    uiwait(fig);
    
    % Callback function for chkBox1
    function chkBox1Changed(chkBox)
        if chkBox.Value
            chkBox2.Enable = 'on';
            chkBox3.Enable = 'on';
        else
            chkBox2.Value = false;
            chkBox3.Value = false;
            chkBox2.Enable = 'off';
            chkBox3.Enable = 'off';
        end
    end
    
    % Callback function for chkBox4
    function chkBox4Changed(chkBox)
        if chkBox.Value
            avgCellDiam.Enable = 'on';
        else
            avgCellDiam.Enable = 'off';
        end
    end
    
    % Function to save values to workspace variables and close GUI
    function saveAndRun()
        % Section 1: Config
        project_dir = longTextBox.Value; % string
        
        % Section 2: Bead Correction
        voxSizeXY = smallTextBox1.Value; % double
        voxSizeZ = smallTextBox2.Value; % double
        
        % Section 2: Checkboxes
        use_bead_buddy = chkBox1.Value; % logical
        doFlatField = chkBox2.Value; % logical
        runAirlocalize = chkBox3.Value; % logical
        
        % Section 3: Image Processing
        make_cell_masks = chkBox4.Value; % logical
        doFlatField_data = chkBox5.Value; % logical
        
        mask_chan = maskChannel.Value; % double
        avg_cell_diameter = avgCellDiam.Value; % double
        
        % Add variables to MATLAB workspace
        assignin('base', 'project_dir', project_dir);
        assignin('base', 'voxSizeXY', voxSizeXY);
        assignin('base', 'voxSizeZ', voxSizeZ);
        assignin('base', 'use_bead_buddy', use_bead_buddy);
        assignin('base', 'doFlatField', doFlatField);
        assignin('base', 'runAirlocalize', runAirlocalize);
        assignin('base', 'make_cell_masks', make_cell_masks);
        assignin('base', 'doFlatField_data', doFlatField_data);
        assignin('base', 'mask_chan', mask_chan);
        assignin('base', 'avg_cell_diameter', avg_cell_diameter);
        
        % Close the GUI
        delete(fig);
        
        % Display confirmation
        disp('Variables saved to workspace. GUI closed.');
        
        % Resume execution of the script
        % Additional actions can be added here if needed
    end
end
