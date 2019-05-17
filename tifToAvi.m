%% One by One Tif to AVI, creates an avi file for each tif file in a directory
% May have very long runtime
% Written by Ethan

function tifToAvi()
    d = uigetdir(pwd, 'Select a folder');
    files = dir(fullfile(d, '*.tif'));
    % Display the names
    A = struct2cell(files);
    
    sizeOfA = size(A);
    numberOfTifs = sizeOfA(1,2);
    
    for k = 1:numberOfTifs
        tifToAviHelper(A{1,k});
%         numberOfFrames = numel(tifInfo);
%         for j = 1:numberOfFrames
%             a = imread(A{1,k}, j);
%             writeVideo(video, a);
%         end
        disp("Done " + k + " / " + numberOfTifs);
    end
end

function tifToAviHelper(tifFile)
    tifInfo = imfinfo(tifFile);
%     tifInfoCells = struct2cell(imfinfo(tifFile));
    justFileName = tifFile(1:end-4);
    name = strcat(justFileName, '.avi');
    video = VideoWriter(name,'Uncompressed AVI');
    
    open(video);
    numberOfFrames = numel(tifInfo);
    for j = 1:numberOfFrames
        a = imread(tifFile, j);
        writeVideo(video, a);
    end
    close(video);
end
