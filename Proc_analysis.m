
% Load in proc file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% Initialize
clearvars
close all

% Session info
dFldrs = {
      'D:\data\M56\20190501_postConditioning\Cond3';...
%          'E:\data\M62\20190511\Cond1';...
%         'E:\data\M62\20190512\Cond2';...
%          'E:\data\M62\20190513\Cond3';...
%         'E:\data\M62\20190515\PostCond2';...
%        'E:\data\M58\20190504\postCond3';...
     
}; % folders containing data
% Other Analsys settings
WL = 1;  %  0/1 = Keep/Remove trials without licks
shkMans = [0 0 0 1 1 0]; % sessions need manual adjustments
shkFlps = {[],[],[],[10 19 28 37 46 55 64 73 82 91],[10 19 28 37 46 55 64 73 82 91],[],[]}; % trials to flip shock on/off
nSesh = length(dFldrs);
% Ryans Figure Settings
RyansFigureSettings
maxLickDiff = [];


%% %% more Initialize
%% Lickport analysis
runSeshs = 1:length(dFldrs); % 5; %
for sesh = runSeshs
    dFldr = dFldrs{sesh};
    seshN = 's';  %Session name
    
   
    shkMan = shkMans(sesh); % used a shock protocol but shocker was turned off for some trials
    shkFlp = shkFlps{sesh}; % list of trials to flip shocks on or off
    
    % Files
    % Trial types
    aTrl = 'adpt';
    sTrl = 'stim';
    % Measurement types
    nCam1 = 'Camera 1';
    nCam2 = 'Camera 2';
    nLick = 'Lick Port';
    nMWhl = 'Mouse Wheel';
    nTrig = 'Trigger';
    nValv = 'Valve';
    nShck = 'Shock';
    nMove = 'Move';
    nProt = 'Protocol';
    nStim = 'angle';
    
    % find files
    dFldr = [dFldr '\'];
    cFNames = num2cell(ls(dFldr),2); %cell array of file names
    fSesh = cellContainsStr(cFNames,seshN); % Session files
    fCam1 = cellContainsStr(cFNames,nCam1); %Camera 1 files
    fCam2 = cellContainsStr(cFNames,nCam2);
    fLick = cellContainsStr(cFNames,nLick);
    fMWhl = cellContainsStr(cFNames,nMWhl);
    fTrig = cellContainsStr(cFNames,nTrig);
    fValv = cellContainsStr(cFNames,nValv);
    fShck = cellContainsStr(cFNames,nShck);
    fMove = cellContainsStr(cFNames,nMove);
    fProt = cellContainsStr(cFNames,nProt);
    fStim = cellContainsStr(cFNames,nStim);
    fTxt = cellContainsStr(cFNames,'txt'); % file type
    fTiff = cellContainsStr(cFNames,'tif');
    seshTitle = cFNames{10}(1:find('_'==cFNames{10},1,'first')-1);

    %% Read in Protocol Settings
    iTOI = fProt;
    fNm = strtrim([dFldr cFNames{iTOI}]);
    protocol = readtable(fNm, 'ReadVariableNames', true, 'Delimiter', '\t');
    
    iTOI = find(fStim);
    fNm = strtrim([dFldr cFNames{iTOI}]);
    stimInfo = load(fNm);
    
    % number of trials
    nTrls = str2num(protocol.TrialName{1});
    nTrlsCheck = length(protocol.TrialName)-1;
    if nTrls~=nTrlsCheck
        warning('Mismatched number of trials in Protocol')
    end
    
    % read trial times (s)
    trlDurs = cellfun(@str2num,protocol.Record(2:end));
    vlvDels = cellfun(@str2num,protocol.FluidValve(2:end));
    vlvDurs = cellfun(@str2num,protocol.Var5(2:end));
    lckDels = cellfun(@str2num,protocol.Lickport(2:end));
    lckDurs = cellfun(@str2num,protocol.Var7(2:end));
    stmDels = stimInfo.staticMovingT(1)*ones(size(trlDurs));
    stmDurs = stimInfo.staticMovingT(2)*ones(size(trlDurs));
    shkDels = cellfun(@str2num,protocol.ShockControl(2:end));
    shkDurs = cellfun(@str2num,protocol.Var10(2:end));        %double check if this is correct
    
    % read from this single trial to get all parameters
    readTrial = 6;
    trlDur = trlDurs(readTrial);
    vlvDel = vlvDels(readTrial);
    vlvDur = vlvDurs(readTrial);
    lckDel = lckDels(readTrial);
    lckDur = lckDurs(readTrial);
    stmDel = stmDels(readTrial);
    stmDur = stmDurs(readTrial);
    
    % Shock times - adjusting for manual changes
    if shkMan==1;
        shkDels(shkFlp) = NaN;
        shkDurs(shkFlp) = NaN;
    end
    shkDels(shkDels==-1) = nan;
    shkDurs(shkDurs==-1) = nan;
    
    shkDel = min(shkDels);
    shkDur = max(shkDurs);
    
    % 0/1 = Don't/Do count spikes that occur during shock time window
    if sum(~isnan(shkDels))>0; shkOn=1; else shkOn=0; end
end

%%  %% Read in Lick data
    
    trls = [1:nTrls];   %trials to process
    lickTimes = cell(size(trls));
    lickCounts = nan(size(trls));
    lickBinCountsT = nan(length(trls),3);
    lickSeq = nan(size(trls))';
    if shkOn==1, shkAdj=shkDur; else shkAdj=0; end %shock adjustment on or off
    
    parfor v = trls %%loop to process trials
        % for v = trls %%loop to process trials
        disp(['Lick Trial_' num2str(v)])
        % select trial
        fSTrl = cellContainsStr(cFNames,[sTrl sprintf('%03d',v)]); %
        
        % Get lick measurement times
        iTOI = find(fSesh & fLick & fSTrl & fTxt);
        if ~isempty(iTOI)
            lickSeq(v) = 1;
            fNm = strtrim([dFldr cFNames{iTOI}]);
            sData = textread(fNm, '%s', 'delimiter', '\n');
            nSamps = length(sData);     % number of samples
            mSecLick = nan(nSamps,1);   % Timestamp of each lick sample
            lLick = nan(nSamps,1);      % Logical of each lick sample
            for u = 1:nSamps
                sStr = sData{u};
                offset = find(sStr=='.',1); % adjusts for change in string length due to more digits
                mSecLick(u) = str2num(sStr(1:offset+6));
                lLick(u) = logical(str2num(sStr(offset+8)));
            end
            sec1Lick = str2num(fNm(end-22:end-15));
            secsLick = mSecLick/1000+ sec1Lick;
            
            % find lick starts
            lickStarts = diff(lLick)==1;
            lickTimes{v} = mSecLick(find(lickStarts)+1)/1000;
            lickCounts(v) = length(lickTimes{v});
        else
            lickSeq(v) = 0;
        end
        % select trial
    end
    
    beep

    %% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %% 

%%%%% Pupil analysis %%%%%%

% Convert vector to cells separated by sessions
nses = length(proc.nframes);
nframes = proc.nframes;
emptycellarray = cell(nses,1);
areavec = proc.pupil.area;

framebegin = 1;
frameend = nframes;

for u = 1:nses;
    emptycell = emptycellarray{u};
    nframe = nframes(u);
    frameend = framebegin + nframe - 1;
    emptycellarray{u} = areavec(framebegin:frameend);
    framebegin = framebegin + nframe;
    areacells = emptycellarray;
    areacells
end
%% 
areamat = areacells;

s = size(areacells);

longest = 0;
for j = 1:s(1)
    m = size(areamat{j});
    if m(1) > longest
        longest = m(1);
    end
end

cellmax = longest;
a = areamat{1};
a = a';

b = cellfun( @(c) [c(:) ; NaN(longest-numel(c),1)], areamat,'un',0);

for k = 2:s(1)
    a = vertcat(a, b{k}');
end
size(a)
a(isnan(a))=0;

areamat = a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Plot pupil area
x = 0:50/1565:16
num_sesplot = 91
plot(x, areacells{num_sesplot})
%xlim([1 500])
xlabel("Time (s)")
ylabel("Pupil area (pixels)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
% for u = 1:nTrls
%     pupilSize = pupilSizes{u};
%     
%     pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
% end

% pupilAbsMax = max(pupilSize500(:));
% % pupilAbsMax = max(abs(pupilSize500(:)));
% pupilAllNorm = pupilSize500/pupilAbsMax;
% % wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
%     compare raw, interpolated and smoothed data
    subplot(rows,cols,1), cla, hold on
    line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
%     plot(wheelTimes{v},wheelCounts{v},'-','color',color)
%     plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
%     plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
%     axis tight, xlim([0 trlDur]),  box off
%     pause
    color = colors{trialType(v)}; 
    pupilTimeFrame = 1:length(pupilSize500(v,:))+1;
    pupilTimeFrame500 = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
%    plot(pupilTimeFrame500,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,areamat(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,areamat(typeTrials,:),'-','color',color)
 %   ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(areamat(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((areamat(typeTrials,:)));
    pupilSte = nanste((areamat(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,areamat(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
% pupilSizes{u} = pupilSize
hold off

%% 
% Generate dt list

pupilareadt = cell(1,nses) ;
for u = 1:nses ; 
    areacell = areacells{u} ;
    t = 1:length(areacell) ;
    v = zeros(length(t)-1,1) ;
    for i = 1:length(t)-1 ;
      v(i) = (areacell(i+1)-areacell(i))/(t(i+1)-t(i)) ;
      pupilareadt{u} = v ;
    end
end
pupilareadt

%% 
% Filter pupil area cell array based on dt

pupilarea_filtered = cell(1,nses);

for index=1:nses
    disp(index)
    area_ses_dt = single(pupilareadt{index}) ;
    area_ses_dt = [NaN;pupilareadt{index}(1:end)];
    
    area_ses = areacells{index};
    %ses1 = table(area_ses1_dt, area_ses1);
    % extract first row: ses1(1:501, 1)
    isgoodframe = (-25 < area_ses_dt & area_ses_dt < 25) ;
    area_ses(~isgoodframe) = nan;
    
    pupilarea_filtered{index} = area_ses;
end

pupilareaProc = pupilarea_filtered;
%% 

pupilarea_filtered_mat = pupilarea_filtered

s = size(pupilarea_filtered);

longest = 0;
for j = 1:s(1)
    m = size(pupilarea_filtered_mat{j});
    if m(1) > longest
        longest = m(1);
    end
end

cellmax = longest;
a = pupilarea_filtered_mat{1};
a = a';

b = cellfun( @(c) [c(:) ; NaN(longest-numel(c),1)], pupilarea_filtered_mat,'un',0);

for k = 2:s(1)
    a = vertcat(a, b{k}');
end
size(a)
%a(isnan(a))=0;

pupilarea_filtered_mat = a

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot original and filtered pupil areas
num_sesplot = 91;   

plot(x, areacells{num_sesplot})
hold on
plot(x, pupilarea_filtered{num_sesplot})
hold off
xlabel("Time (s)")
ylabel("Pupil area (pixels)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
% for u = 1:nTrls
%     pupilSize = pupilSizes{u};
%     
%     pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
% end

% pupilAbsMax = max(pupilSize500(:));
% % pupilAbsMax = max(abs(pupilSize500(:)));
% pupilAllNorm = pupilSize500/pupilAbsMax;
% % wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
%     compare raw, interpolated and smoothed data
    subplot(rows,cols,1), cla, hold on
    line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
%     plot(wheelTimes{v},wheelCounts{v},'-','color',color)
%     plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
%     plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
%     axis tight, xlim([0 trlDur]),  box off
%     pause
    color = colors{trialType(v)}; 
    pupilTimeFrame = 1:length(pupilSize500(v,:))+1;
    pupilTimeFrame500 = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
%    plot(pupilTimeFrame500,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,pupilarea_filtered_mat(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,pupilarea_filtered_mat(typeTrials,:),'-','color',color)
 %   ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(pupilarea_filtered_mat(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((pupilarea_filtered_mat(typeTrials,:)));
    pupilSte = nanste((pupilarea_filtered_mat(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,pupilarea_filtered_mat(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
% pupilSizes{u} = pupilSize
hold off
%% 

    com = proc.pupil.com;
    
    area_raw = proc.pupil.area_raw;
    area = proc.pupil.area;
    
    figure;
    
    x0 = round(proc.locROI{1,5}(1) * proc.sc) ;
    y0 = round(proc.locROI{1,5}(2) * proc.sc);
    x1 = round(x0 + proc.locROI{1,5}(3) * proc.sc);
    y1 = round(y0 +  proc.locROI{1,5}(4) * proc.sc);
    disp("x0: " + x0)
    disp("y0: " + y0)
    disp("x1: " + x1)
    disp("y1: " + y1)
    
    
    
    for k = 1:proc.nframes(1,1)
        
        title("Frame:" + k);
        a = imread("Cond3_00000_stim001_05_01_2019_15_30_43.30797__Camera 1.tif", k);
        b=a(y0:y1,x0:x1);
        OutputImage = imadjust(b,[0; 1],[.4; 1])
        
        subplot(1,4,1), imshow(OutputImage);
        title("video");
        
        subplot(1, 4, 2), viscircles([0 0], sqrt(area_raw(k)/pi));
        title("raw data");
        
        subplot(1, 4, 3), viscircles([0 0], sqrt(area(k)/pi));
        title("processed data");
        
        subplot(1, 4, 4), viscircles([0 0], sqrt(pupilarea_filtered{1,1}(k)/pi));
        title("processed filtered   data");
        pause(0.1);
        clf;
    end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% CoM Analysis -- X component %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert vector to cells separated by sessions

emptycellarray = cell(nses,1);
comvec = proc.pupil.com;
comvecX = comvec(1:end, 1)
comvecY = comvec(1:end, 2)

framebegin = 1;
frameend = nframes;

for u = 1:nses;
    emptycell = emptycellarray{u};
    nframe = nframes(u);
    frameend = framebegin + nframe - 1;
    emptycellarray{u} = comvecX(framebegin:frameend);
    framebegin = framebegin + nframe;
    comcellsX = emptycellarray;
    comcellsX
end
%% 
comxmat = comcellsX

s = size(comcellsX);

longest = 0;
for j = 1:s(1)
    m = size(comxmat{j});
    if m(1) > longest
        longest = m(1);
    end
end

cellmax = longest;
a = comxmat{1};
a = a';

b = cellfun( @(c) [c(:) ; NaN(longest-numel(c),1)], comxmat,'un',0);

for k = 2:s(1)
    a = vertcat(a, b{k}');
end
size(a)
a(isnan(a))=0;

comxmat = a

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Plot pupil area

num_sesplot = 91
plot(x, comcellsX{num_sesplot})
ylim([82 86])
xlim([1 16])
xlabel("Time (s)")
ylabel("Pupil position (pixels)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
% for u = 1:nTrls
%     pupilSize = pupilSizes{u};
%     
%     pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
% end

% pupilAbsMax = max(pupilSize500(:));
% % pupilAbsMax = max(abs(pupilSize500(:)));
% pupilAllNorm = pupilSize500/pupilAbsMax;
% % wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
%     compare raw, interpolated and smoothed data
    subplot(rows,cols,1), cla, hold on
    line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
%     plot(wheelTimes{v},wheelCounts{v},'-','color',color)
%     plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
%     plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
%     axis tight, xlim([0 trlDur]),  box off
%     pause
    color = colors{trialType(v)}; 
    pupilTimeFrame = 1:length(pupilSize500(v,:))+1;
    pupilTimeFrame500 = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
%    plot(pupilTimeFrame500,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,comxmat(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,comxmat(typeTrials,:),'-','color',color)
 %   ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(comxmat(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((comxmat(typeTrials,:)));
    pupilSte = nanste((comxmat(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,comxmat(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
% pupilSizes{u} = pupilSize
hold off


%% 
% Generate dt list

pupilcomdtX = cell(1,nses) ;
for u = 1:nses ; 
    comcellX = comcellsX{u} ;
    t = 1:length(comcellX) ;
    v = zeros(length(t)-1,1) ;
    for i = 1:length(t)-1 ;
      v(i) = (comcellX(i+1)-comcellX(i))/(t(i+1)-t(i)) ;
      pupilcomdtX{u} = v ;
    end
end
pupilcomdtX

%% 
% Filter pupil movement cell array based on dt

pupilcomX_filtered = cell(1,nses);

for index=1:nses
    disp(index)
    com_sesX_dt = single(pupilcomdtX{index}) ;
    com_sesX_dt = [NaN;pupilcomdtX{index}(1:end)];
    
    com_ses = comcellsX{index};
    %ses1 = table(area_ses1_dt, area_ses1);
    % extract first row: ses1(1:501, 1)
    isgoodframe = (-1 < com_sesX_dt & com_sesX_dt < 1) ;
    com_ses(~isgoodframe) = nan;
    
    pupilcomX_filtered{index} = com_ses;
end

pupilcomProc = pupilcomX_filtered;    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot original and filtered pupil movement

num_sesplot = 3
plot(comcellsX{num_sesplot})
hold on
plot(pupilcomX_filtered{num_sesplot})
hold off

%% 

pupilcomX_filtered_mat = pupilcomX_filtered

s = size(comcellsX);

longest = 0;
for j = 1:s(1)
    m = size(pupilcomX_filtered_mat{j});
    if m(1) > longest
        longest = m(1);
    end
end

cellmax = longest;
a = pupilcomX_filtered_mat{1};
a = a';

b = cellfun( @(c) [c(:) ; NaN(longest-numel(c),1)], pupilcomX_filtered_mat,'un',0);

for k = 2:s(1)
    a = vertcat(a, b{k}');
end
size(a)
% a(isnan(a))=0;

pupilcomX_filtered_mat = a


%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
% for u = 1:nTrls
%     pupilSize = pupilSizes{u};
%     
%     pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
% end

% pupilAbsMax = max(pupilSize500(:));
% % pupilAbsMax = max(abs(pupilSize500(:)));
% pupilAllNorm = pupilSize500/pupilAbsMax;
% % wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
%     compare raw, interpolated and smoothed data
    subplot(rows,cols,1), cla, hold on
    line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
%     plot(wheelTimes{v},wheelCounts{v},'-','color',color)
%     plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
%     plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
%     axis tight, xlim([0 trlDur]),  box off
%     pause
    color = colors{trialType(v)}; 
    pupilTimeFrame = 1:length(pupilSize500(v,:))+1;
    pupilTimeFrame500 = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
%    plot(pupilTimeFrame500,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,pupilcomX_filtered_mat(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,pupilcomX_filtered_mat(typeTrials,:),'-','color',color)
 %   ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(pupilcomX_filtered_mat(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((pupilcomX_filtered_mat(typeTrials,:)));
    pupilSte = nanste((pupilcomX_filtered_mat(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,pupilcomX_filtered_mat(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
% pupilSizes{u} = pupilSize
hold off

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% CoM Analysis -- Y component %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert vector to cells separated by sessions

emptycellarray = cell(nses,1);
comvec = proc.pupil.com;
comvecX = comvec(1:end, 1)
comvecY = comvec(1:end, 2)

framebegin = 1;
frameend = nframes;

for u = 1:nses;
    emptycell = emptycellarray{u};
    nframe = nframes(u);
    frameend = framebegin + nframe - 1;
    emptycellarray{u} = comvecY(framebegin:frameend);
    framebegin = framebegin + nframe;
    comcellsY = emptycellarray;
    comcellsY
end

%% 
comymat = comcellsY

s = size(comcellsY);

longest = 0;
for j = 1:s(1)
    m = size(comymat{j});
    if m(1) > longest
        longest = m(1);
    end
end

cellmax = longest;
a = comymat{1};
a = a';

b = cellfun( @(c) [c(:) ; NaN(longest-numel(c),1)], comymat,'un',0);

for k = 2:s(1)
    a = vertcat(a, b{k}');
end
size(a)
a(isnan(a))=0;

comymat = a

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Plot pupil area
num_sesplot = 91
plot(x, comcellsY{num_sesplot})
xlim([1 16])
xlabel("Time (s)")
ylabel("Pupil position (pixels)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
% for u = 1:nTrls
%     pupilSize = pupilSizes{u};
%     
%     pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
% end

% pupilAbsMax = max(pupilSize500(:));
% % pupilAbsMax = max(abs(pupilSize500(:)));
% pupilAllNorm = pupilSize500/pupilAbsMax;
% % wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
%     compare raw, interpolated and smoothed data
    subplot(rows,cols,1), cla, hold on
    line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
%     plot(wheelTimes{v},wheelCounts{v},'-','color',color)
%     plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
%     plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
%     axis tight, xlim([0 trlDur]),  box off
%     pause
    color = colors{trialType(v)}; 
    pupilTimeFrame = 1:length(pupilSize500(v,:))+1;
    pupilTimeFrame500 = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
%    plot(pupilTimeFrame500,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,comymat(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,comymat(typeTrials,:),'-','color',color)
 %   ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(comymat(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((comymat(typeTrials,:)));
    pupilSte = nanste((comymat(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,comymat(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
% pupilSizes{u} = pupilSize
hold off

%% 
% Generate dt list

pupilcomdtY = cell(1,nses) ;
for u = 1:nses ; 
    comcellY = comcellsY{u} ;
    t = 1:length(comcellY) ;
    v = zeros(length(t)-1,1) ;
    for i = 1:length(t)-1 ;
      v(i) = (comcellY(i+1)-comcellY(i))/(t(i+1)-t(i)) ;
      pupilcomdtY{u} = v ;
    end
end
pupilcomdtY

%% 
% Filter pupil area cell array based on dt

pupilcomY_filtered = cell(1,nses);

for index=1:nses
    disp(index)
    com_sesY_dt = single(pupilcomdtY{index}) ;
    com_sesY_dt = [NaN;pupilcomdtY{index}(1:end)];
    
    com_ses = comcellsY{index};
    %ses1 = table(area_ses1_dt, area_ses1);
    % extract first row: ses1(1:501, 1)
    isgoodframe = (-.5 < com_sesY_dt & com_sesY_dt < .5) ;
    com_ses(~isgoodframe) = nan;
    
    pupilcomY_filtered{index} = com_ses;
end
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot original and filtered pupil movement

num_sesplot = 3
plot(comcellsY{num_sesplot})
hold on
plot(pupilcomY_filtered{num_sesplot})
hold off

%% 
pupilcomY_filtered_mat = pupilcomY_filtered

s = size(comcellsX);

longest = 0;
for j = 1:s(1)
    m = size(pupilcomY_filtered_mat{j});
    if m(1) > longest
        longest = m(1);
    end
end

cellmax = longest;
a = pupilcomY_filtered_mat{1};
a = a';

b = cellfun( @(c) [c(:) ; NaN(longest-numel(c),1)], pupilcomY_filtered_mat,'un',0);

for k = 2:s(1)
    a = vertcat(a, b{k}');
end
size(a)
% a(isnan(a))=0;

pupilcomY_filtered_mat = a


%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
% for u = 1:nTrls
%     pupilSize = pupilSizes{u};
%     
%     pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
% end

% pupilAbsMax = max(pupilSize500(:));
% % pupilAbsMax = max(abs(pupilSize500(:)));
% pupilAllNorm = pupilSize500/pupilAbsMax;
% % wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
%     compare raw, interpolated and smoothed data
    subplot(rows,cols,1), cla, hold on
    line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
%     plot(wheelTimes{v},wheelCounts{v},'-','color',color)
%     plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
%     plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
%     axis tight, xlim([0 trlDur]),  box off
%     pause
    color = colors{trialType(v)}; 
    pupilTimeFrame = 1:length(pupilSize500(v,:))+1;
    pupilTimeFrame500 = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
%    plot(pupilTimeFrame500,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,pupilcomY_filtered_mat(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,pupilcomY_filtered_mat(typeTrials,:),'-','color',color)
 %   ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(pupilcomY_filtered_mat(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((pupilcomY_filtered_mat(typeTrials,:)));
    pupilSte = nanste((pupilcomY_filtered_mat(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,pupilcomY_filtered_mat(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
% pupilSizes{u} = pupilSize
hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERPOLATE area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to row vector
A = pupilarea_filtered_mat;
B = reshape(A.',1,[])

% Interp
%// Indices of NaNs

t1 = find(isnan(factor)); 

%// Indices of non-NaNs
t2 = find(~isnan(factor));

%// Get index for each NaN index that is closest, with a tie-case 
%// (closest non-NaN number being at equal distance on either side) 
%// selecting the left one
[~,ind1] = min(abs(bsxfun(@minus,t1,t2'))); %//'

%// Replace NaNs with the closest non-NaNs
factor(t1) = factor(t2(ind1))
%% 

% Outer for loop to iterate between different videos
for u = 1:length(pupilarea_filtered)
    
   
    %Inner for loop to iterate between different frames seraching for a Nan
    for k = 1:length(pupilarea_filtered{1, u})
        %counters to increment up or down frames to find the nearrest
        %values in the while loops
        j = 0;
        m = 0;
        %Variables to save the nearest values in
        closest_prev = 0;
        closest_next = 0;
        
        if isnan(pupilarea_filtered{1, u}(k))
            
            %find previous frame area 
            while k - j > 0
                if (k - j) == 1
                    closest_prev = NaN;
                    break;
                end
                if ~isnan(pupilarea_filtered{1, u}(k - j))
                    closest_prev = pupilarea_filtered{1, u}(k-j);
                end
                j = j + 1;
            end
            
            %find next frame area
            while k + m < (length(pupilarea_filtered{1, u}) + 1)
                if (k + m) == length(pupilarea_filtered{1, u})
                    closest_next = NaN;
                    break;
                end
                if ~isnan(pupilarea_filtered{1, u}(k + m))
                    closest_next = pupilarea_filtered{1, u}(k+m);
                    break;
                end
                m = m + 1;
            end
            
            %average closest_next and closest_prev to get value for
            %previous nan valued frame
            if isnan(closest_prev)
                average = closest_next;
            elseif isnan(closesy_next)
                average = closest_prev;
            else
                average = closest_prev + closest_next;
            end
            pupilarea_filtered{1, u}(k) = average;
        end
    end
   
end


    
    %% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot original and filtered pupil areas
num_sesplot = 4;   

hold off
plot(areacells{num_sesplot})
hold on  
plot(areacells_interp{num_sesplot} - 300)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pupil plots original 

figure(7), clf
figure(6), clf

rows = 1;
cols = 4;

%stim sequence
    stimSeqRaw = stimInfo.contrastSeq;
    if isfield(stimInfo,'wavSeq')
        if sum(stimInfo.wavSeq~=1)>0
            stimSeqRaw = stimInfo.wavSeq;
        end
    end

%from autolog
stimSeq0 = zeros(nTrls,1);
stimSeq0(stimSeqRaw>1) = stimInfo.angleSeq(stimSeqRaw>1);
[~,~,stimId] = unique(stimSeq0);
stimSeq = stimId-1;
% stimSeq = zeros(nTrls,1);
% Trial Types
trialType = nan(size(trls));
trialType(lckDels==-1 & ~stimSeq) = 1; % No Lickport, Blank
trialType(lckDels==-1 & stimSeq) = 2; % No Lickport, Stimulus
trialType(lckDels~=-1 & ~stimSeq) = 3; % Lickport, Blank
trialType(lckDels~=-1 & stimSeq) = 4; % Lickport, Stimulus
nTrialTypes = 4;
%
colors = {'k' 'g' 'b' 'c'};


pupilSize500 = nan(nTrls,500);
% for u = 1:nTrls
%     pupilSize = pupilSizes{u};
%     
%     pupilSize500(u,1:length(pupilSize)) = (pupilSize-pupilSize(1))./pupilSize;
% end

% pupilAbsMax = max(pupilSize500(:));
% % pupilAbsMax = max(abs(pupilSize500(:)));
% pupilAllNorm = pupilSize500/pupilAbsMax;
% % wheelFreqTrlNorm = wheelFreqSmth./repmat(max(abs(wheelFreqSmth)')',1,length(wheelFreqSmth));

trlSclr = .1;
for v = trls
%     compare raw, interpolated and smoothed data
    subplot(rows,cols,1), cla, hold on
    line([0 trlDur],[0 0],'color',[.5 .5 .5],'linestyle','--')
%     plot(wheelTimes{v},wheelCounts{v},'-','color',color)
%     plot(wheelTimeMs,wheelCountMs(v,:),'--','color','r')
%     plot(wheelTimeMs,wheelCountSmth(v,:),'--','color','g')
%     axis tight, xlim([0 trlDur]),  box off
%     pause
    color = colors{trialType(v)}; 
    pupilTimeFrame = 1:length(pupilSize500(v,:))+1;
    pupilTimeFrame500 = 1:length(pupilSize500(v,:));

    subplot(5,4,4*4+[1 2]), hold all
%    plot(pupilTimeFrame500,pupilSize500(v,:),'-','color',color)
    axis tight
    
    subplot(1,4,3), hold all
    plot(pupilTimeFrame,areamat(v,:)+(v-1)*trlSclr,'-','color',color)
    axis tight
end


figure(7), clf
legendTypes = {'No Lickport, Blank'; 'No Lickport, Stimulus'; 'Lickport, Blank'; 'Lickport, Stimulus'};
c=0;
for u = 1:nTrialTypes
    color = colors{u};
    % Grouped plots
    figure(6), subplot(5,nTrialTypes,nTrialTypes*(u-1)+[1 2]), hold on
    typeTrials = find(trialType==u);
    plot(pupilTimeFrame,areamat(typeTrials,:),'-','color',color)
 %   ylim([-1 1])
    legend(legendTypes{u},'location','southeast')
    
    
    figure(7), hold on
    % Imagesc plots
    subplot(nTrialTypes,2,(u-1)*2+1)
    imagesc('xdata',pupilTimeFrame,'cdata',abs(areamat(typeTrials,:)))
    colormap parula
    caxis([0 1]), axis tight
    title(legendTypes{u})
    % PSTH Overlays
    subplot(nTrialTypes,2,[2:2:8]), hold on
    pupilMean = nanmean((areamat(typeTrials,:)));
    pupilSte = nanste((areamat(typeTrials,:)));
%     pupilSte(isnan(pupilSte)) = mean(isnan(pupilSte)-1)
    patch([pupilTimeFrame fliplr(pupilTimeFrame)], [pupilMean+pupilSte fliplr(pupilMean-pupilSte)],color,'FaceAlpha',.1,'edgecolor','none')
    plot(pupilTimeFrame,pupilMean,'-','color',color)
    
    
    figure(6), subplot(1,nTrialTypes,nTrialTypes), hold all
    for v = typeTrials
        plot(pupilTimeFrame,areamat(v,:)+c,'-','color',color);
        c=c+.1;
    end
    c=c+1;
end
% pupilSizes{u} = pupilSize
hold off
%% 

   com = proc.pupil.com;
    
    area_raw = proc.pupil.area_raw;
    area = proc.pupil.area;
    
    figure;
    hold on
    x0 = round(proc.locROI{1,5}(1) * proc.sc) ;
    y0 = round(proc.locROI{1,5}(2) * proc.sc);
    x1 = round(x0 + proc.locROI{1,5}(3) * proc.sc);
    y1 = round(y0 +  proc.locROI{1,5}(4) * proc.sc);
    disp("x0: " + x0)
    disp("y0: " + y0)
    disp("x1: " + x1)
    disp("y1: " + y1)
    

    
    for k = 1:proc.nframes(1,1)
        
        
        
        
        a = imread("Cond3_00000_stim001_05_01_2019_15_30_43.30797__Camera 1.tif", k);
        b=a(y0:y1,x0:x1);
        
        
        
%         subplot(1,4,1), imshow(b, 'InitialMagnification', 100);
%         title("Video");
        
        subplot(1, 3, 1),
        
        I = image(b);
        hold on;
        circle(com(k, 1),com(k, 2),sqrt(area_raw(k)/pi));
        title("Raw Data");
        
        
        subplot(1, 3, 2), 
        
        J = image(b);
        colormap(gray(256));
        hold on;
        circle(com(k,1),com(k,2),sqrt(area(k)/pi));
        title("Processed Data");
        
        subplot(1, 3, 3), 
        
        K = image(b);
        hold on;
        circle(com(k,1),com(k,2),sqrt(pupilarea_filtered{1,1}(k)/pi));
        title("Processed In terpolated Data");
        sgtitle("Pupil Area Plots, Frame: " + k);
        pause(0.2);
        clf;
    end
% Circle function better than viscircles
function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'r', 'LineWidth', 2);
    uistack(h, 'top');
    hold off
end