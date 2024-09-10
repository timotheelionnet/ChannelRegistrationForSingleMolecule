

function [] = verify_channel_key(key_tab)
% make sure only one ch desiganted as ref
if sum(key_tab.isReference) ~= 1
    disp('Please revise your channel key')
    error('You must designate only one channel as your reference channel')
    
end

% make sure ref channel is localizable
if key_tab.isLocalizable( find(key_tab.isReference)) ~= 1 
    disp('Please revise your channel key')
    error('Your reference channel must be localizable')
end
end
