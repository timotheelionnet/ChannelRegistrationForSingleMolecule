function imSize = getImageSizeFromFile(fileName)
    info = imfinfo(fileName);
    imSize(1) = info(1).Width;
    imSize(2) = info(1).Height;
    z = numel(info);
    if z >1
        imSize(3) = z;
    end
end