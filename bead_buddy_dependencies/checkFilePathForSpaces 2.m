

function checkFilePathForSpaces(filePath)
    % checkFilePathForSpaces Checks if the specified file path contains spaces.
    %   If the file path contains spaces, it returns an error message.
    %   If the file path does not contain spaces, it confirms the path is valid.

    if contains(filePath, ' ')
        disp(filePath)
        error('Error: This filepath cannot contain spaces');

    end
end