function plotPCs(h, pc)
    R = cast(h.nframes, 'like', 5);
    N = pc;
    X = linspace(1, R, R);
    Y1 = h.motSVD{1,1}(:,N);
    
    % subplot(m,n,p) divides the current figure into an m-by-n grid and
    % creates axes in the position specified by p. MATLAB® numbers subplot
    % positions by row. The first subplot is the first column of the first
    % row, the second subplot is the second column of the first row, and
    % so on. If axes exist in the specified position, then this command 
    % makes the axes the current axes.
    
    figure;
    subplot(2,1,1);
    plot(X,Y1)
    title("PC #" + pc)
    
    implay(h.files{1,1}, 5);
   
end