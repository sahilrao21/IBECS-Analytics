%% Description of what this function does:
% This function concatenates any number of tif files, saves the concatenated 
% tif file in the current folder, and converts and saves the concatenated 
% tif file as an AVI file for use in MovieGUI.
%
% Possible use: to combine short trials to use in MovieGUI
%
% Arguments: 
%   tifFileArray: a cell array of tif files that you want to concatenate 
%                 and convert into avi (ex. {'tif1.tif', 'tif2.tif'})
%   destinationTif: what you want the name of the concatenated tif file to
%                   be (ex. 'destination.tif'). It will be saved in the
%                   same folder as this script is.
%   destinationAVI: what you want the name of the concatenated tiff
%                   converted into an AVI to be. (ex. 'destination.avi'). 
%                   It will be saved in the same folder as this script.
%
% Written by Ethan Mehta

function tifConcatenate(tifFileArray, destinationTif, destinationAVI)
    total_images = 0;
    for k = 1:length(tifFileArray)
        total_images = total_images + numel(imfinfo(tifFileArray{k}));
    end
    disp(total_images)
    
    A = cell(1, total_images);
    
    disp(length(tifFileArray));
    
    image_counter = 1;
    for k = 1:length(tifFileArray)
        num_images = numel(imfinfo(tifFileArray{k}));
        
        for j = image_counter:image_counter + num_images - 1
            A{1,j} = imread(tifFileArray{k}, j - image_counter + 1);
        end
        
        image_counter = image_counter + num_images;
    end
    
    imwrite(A{1, 1}, destinationTif);
    for k = 2:total_images
        imwrite(A{1, k}, destinationTif,'WriteMode','append');
    end
 
    % Convert concatenated tif to an avi
     v = VideoWriter(destinationAVI,'Uncompressed AVI');
     open(v);
     for k=1:total_images      % assumes 10 images to write to file
         writeVideo(v, A{1,k});
     end
     close(v);

    % Remove Transition Frames
    % In Progress if needed
end

