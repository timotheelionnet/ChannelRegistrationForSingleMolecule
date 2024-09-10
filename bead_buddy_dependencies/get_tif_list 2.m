


% retruns list of paths to all tifs in input dir
function img_path_list = get_tif_list(my_dir)

all_files = dir(fullfile(my_dir, '*.tif'));

imgs_tab = struct2table(all_files);

img_paths = strings(size(imgs_tab.folder));

for i =1:numel(img_paths)
    
    cur_path = fullfile( string(imgs_tab.folder(i)), string(imgs_tab.name(i)) );
    
    img_paths(i) = string(cur_path);
end
    
img_path_list = img_paths;

end
