% What this function does: Let's you crop the *.tif files to only include
% polygon sections that you include in a drawing. It is much easier to use
% this than using the select ROI tools provided in MovieGUI. 
% Arguments: tifFile: the name of the *.tif file that you want to select a
%                     region of.
%
%            destinationTif: the name of the file containing the cropped
%                            image file.
% After running this function use the cropped image *.tif file in
% MovieGUI.

function selectROIimproved(tifFile, destinationTif) 
    % Display an image
    figure;
    imshow(imread(tifFile, 1));
    a = imfinfo(tifFile);

    % Begin interactive placement of a polygon
    h = drawpolygon();

    % Make the polygon face transparent but still clickable
    h.FaceAlpha = 0;
    h.FaceSelectable = true;

    x = round(h.Position(1:0.5*end));
    y = round(h.Position(0.5*end+1:end));
    
    minx = max(min(x), 1);
    maxx = min(max(x), 1024);

    miny = max(min(y), 1);
    maxy = min(max(y), 768);
    
    bw = poly2mask(x,y,768,1024);
    
    % Uncomment this section to see the mask defined by the polygon section.
    % imshow(bw);
    % hold on
    % plot(x,y,'b','LineWidth',2)
    % hold off
    % imshow(bw);
    
    for k = 1:numel(a)
        mask_applied = uint8(bw) .* imread(tifFile, k);
        imwrite(mask_applied(miny:maxy, minx:maxx), destinationTif,'WriteMode','append');
        disp(k)
    end
end











