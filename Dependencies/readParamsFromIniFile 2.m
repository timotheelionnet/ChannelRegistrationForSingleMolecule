function params = readParamsFromIniFile(configFileName)
% reads all parameters (except file patterns in [Files] section from
% iniconfig file configFileName.
% parameters are saved in structure params with the following pattern:
%params.<section name>.<key name>


% making sure the ini file exists
if ~exist(configFileName,'file')
    params = [];
    disp(['Could not find the ini config file: ',configFileName]);
    return
end

% load ini config file
inConf = IniConfig();
inConf.ReadFile(configFileName);
sections = inConf.GetSections();

% populate parameter structure with all entries outside of [Files] section
params = struct();
for i=1:numel(sections)
    if ~strcmp(sections{i},'[Files]')
        sections{i} = strrep(sections{i},'[','');
        sections{i} = strrep(sections{i},']','');
        
        if ~isfield(params,sections{i})
            params.(sections{i}) = [];
        end
        
        [keys, count_keys] = inConf.GetKeys(sections{i});
        
        for j=1:numel(keys)
            params.(sections{i}).(keys{j}) = inConf.GetValues(sections{i}, keys{j});
        end
    end
    
end





