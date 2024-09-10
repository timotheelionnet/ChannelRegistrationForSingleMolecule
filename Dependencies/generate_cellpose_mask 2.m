function mask = generate_cellpose_mask(image_path, avg_cell_diameter, my_thresh, save_dir)


% lets try using cellpose from matlab with the medical imgaging toolbox

imgpath = image_path;

disp('Processing');
disp(imgpath);

img = imread(imgpath);

if ~exist(save_dir, "dir")
    mkdir(save_dir);
    disp('Created dir for saving')
    disp(save_dir)
end

% show image
% figure
% imshow(img)

img_bw = im2gray(img);

%% show in grays (if not already there)
% imshow(img_bw)


averageCellDiameter = avg_cell_diameter; 

% make cellpos model, depends on Medical Imaging Toolbox Cellpose
% Interface, Deep Learning Toolbox, Computer Vision Toolbox, Parallel
% computing toolbox (for use of GPU)
cp = cellpose(Model="nuclei", ExecutionEnvironment = "gpu", Acceleration="none");

%%
labels = segmentCells2D(cp,img_bw,ImageCellDiameter=averageCellDiameter, CellThreshold= my_thresh);


%%
figure
tiledlayout(1,2,TileSpacing="none",Padding="tight")
nexttile
imshow(img, Colormap = gray, DisplayRange = [50, quantile(img(:), 0.98)])
title("raw Image")
nexttile
imshow(labeloverlay(img,labels))
title("Predicted Labels")

% linkaxes(findobj(gcf,Type="axes"));

% close all
%%
mask = labels +1; % necsesary to add 1 because tif writing tries to index usto 0

[~, n] = fileparts(imgpath);

disp('found cells')
disp(numel(unique(labels)) -1)

%%

save_filepath = fullfile(save_dir, strcat(n, '_MASK.tif')); 

imwrite(mask, gray, save_filepath, 'tiff');

disp('Mask saved');
disp(save_filepath);
disp('~~~~~~~~~~~~~')

end