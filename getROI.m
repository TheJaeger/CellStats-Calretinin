function ROI = getROI(imageArray)
%   Author: Siddhartha Dhiman
%   e-mail: sdhiman@buffalo.edu
%   -----------------------------------------------------------------------
%   getROI.m computes circular ROIs and stores them in the 3D
%   matrix 'ROI' defined by MxNxP. MxN represents the image size and P
%   represents the number of ROIs discovered.
%   -----------------------------------------------------------------------
%   Input Arguments
%       imageArray: Image file annotated with circualr ROI
%   -----------------------------------------------------------------------
%   Output
%       ROI: Binary image mask with ROI. Multiply this by image to extract
%       contents within the ROI.
%   -----------------------------------------------------------------------

%%   Extract blue channel
IB = imbinarize(imageArray(:,:,3));
IB_fill = imfill(IB, 'holes');

%% Locate Objects
regsIB = regionprops(IB_fill, 'Centroid', 'Area');

%% Find and Remove Smallest Element (Line)
[val,ind] = min([regsIB.Area]);
regsIB(ind) = [];

%% Find Radius of Each ROI
for r = 1:length(regsIB)
    R(r) = sqrt(regsIB(r).Area/pi);
end

%% Reform ROIs
imageSizeX = size(IB,2);
imageSizeY = size(IB,1);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
for s = 1:size(regsIB, 1)
    ROI(:,:,s) = (rowsInImage - regsIB(s).Centroid(2)).^2 ... 
        + (columnsInImage - regsIB(s).Centroid(1)).^2 <= R(s).^2;
end

clearvars variables -except ROI


