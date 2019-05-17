%% Description of what this function does:
% Directs you to choose a directory, picks out all of the tif files in the
% directory, and makes them into one large tif file 1 by 1,
% non-destructively. 
% Written by Ethan Mehta

function tifToAviConcatenate()
    d = uigetdir(pwd, 'Select a folder');
    files = dir(fullfile(d, '*.tif'));
    % Display the names
    A = struct2cell(files);
    
    sizeOfA = size(A);
    numberOfTifs = sizeOfA(1,2);
    
    video = VideoWriter('concat1.avi','Uncompressed AVI');
    open(video);
    for k = 1:numberOfTifs
        tifInfo = imfinfo(A{1,k});
        numberOfFrames = numel(tifInfo);
        for j = 1:numberOfFrames
            a = imread(A{1,k}, j);
            writeVideo(video, a);
        end
        disp("Done " + k + " / " + numberOfTifs);
    end
    close(video);
end


