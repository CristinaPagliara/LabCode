clear
clc

% 1 import and set
expRootFolder = 'E:\PhD\E04_PERIPE_PP-CP\Data\ephys analyses';
load(fullfile(expRootFolder, 'pooledTable.mat'), 'T')
load(fullfile(expRootFolder, 'pooledDataSet.mat'), 'ts')


%% select sub-tables
for idxMouse = 1:size(ts,2)
    eventTimes(idxMouse).mouse = ts(idxMouse).T(:, ["sideIdx", "endloc", "speeds", "phase", "isDarkTrial", "Ton", "Toff", "Ton_rec", "Toff_rec"]);
    spikeTimes(idxMouse).mouse = T(T.mouseNum == idxMouse, ["kID", "acronym", "templWaveform", "spikeTimes"]);
end
clear T ts


%% Extract units'spiketimes in a given region of a given mouse
% idxMouse = 8;
% idxArea = 'CA1';
% mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
% st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
% 
% 
% 
% timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
% stimuli = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.phase == 2));
% 
% for idxUnit = 1:size(spikeTimes,2)
%     spikes = st{idxUnit};
%     [raster,indices] = Sync(spikes, stimuli, 'durations', [-1 1]);     % compute spike raster data
%     [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
%     figure;
%     tiledlayout(2,1);
%     nexttile
%     PlotSync(raster,indices, 'durations', [-1 1]);
%     nexttile
%     plot(t, m);
% end

%% striatum slow contra
idxArea = 'CP';
psthPRE = [];
psthPOST = [];
auROC_PRE = [];
auROC_POST = [];
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPRE = [psthPRE m];
            [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPOST = [psthPOST m];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_PRE)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_PRE = [auROC_PRE; auROC];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_POST)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_POST = [auROC_POST; auROC];

        end
    end
end

zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
xlabel('time (s)')
ylabel('zscored firing rate')
[h p] = ttest(auROC_PRE, auROC_POST)
nexttile
plot(auROC_PRE, auROC_POST, 'o')
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'CP slow', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
setCurrentFigureSize(937, 420)


%% striatum fast contra
idxArea = 'CP';
psthPRE = [];
psthPOST = [];
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPRE = [psthPRE m];
            [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPOST = [psthPOST m];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_PRE)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_PRE = [auROC_PRE; auROC];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_POST)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_POST = [auROC_POST; auROC];

        end
    end
end

zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
xlabel('time (s)')
ylabel('zscored firing rate')
[h p] = ttest(auROC_PRE, auROC_POST)
nexttile
plot(auROC_PRE, auROC_POST, 'o')
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'CP fast', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
setCurrentFigureSize(937, 420)

%% CA1 slow contra
idxArea = 'CA1';
psthPRE = [];
psthPOST = [];
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPRE = [psthPRE m];
            [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPOST = [psthPOST m];
            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_PRE)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_PRE = [auROC_PRE; auROC];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_POST)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_POST = [auROC_POST; auROC];
        end
    end
end

zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
xlabel('time (s)')
ylabel('zscored firing rate')
[h p] = ttest(auROC_PRE, auROC_POST)
xlabel('time (s)')
ylabel('zscored firing rate')
nexttile
plot(auROC_PRE, auROC_POST, 'o')
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'CA1 slow', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
setCurrentFigureSize(937, 420)
%% CA1 fast contra
idxArea = 'CA1';
psthPRE = [];
psthPOST = [];
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPRE = [psthPRE m];
            [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPOST = [psthPOST m];
            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_PRE)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_PRE = [auROC_PRE; auROC];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_POST)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_POST = [auROC_POST; auROC];
        end
    end
end

zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
[h p] = ttest(auROC_PRE, auROC_POST)
xlabel('time (s)')
ylabel('zscored firing rate')
nexttile
plot(auROC_PRE, auROC_POST, 'o')
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'CA1 fast', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
setCurrentFigureSize(937, 420)

%% CSSp-ul6a slow contra
idxArea = 'SSp-ul6a';
psthPRE = [];
psthPOST = [];
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 1) & (timestamps.isDarkTrial == 0));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPRE = [psthPRE m];
            [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPOST = [psthPOST m];
            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_PRE)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_PRE = [auROC_PRE; auROC];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_POST)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_POST = [auROC_POST; auROC];
        end
    end
end

zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
[h p] = ttest(auROC_PRE, auROC_POST)
xlabel('time (s)')
ylabel('zscored firing rate')
nexttile
plot(auROC_PRE, auROC_POST, 'o')
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'SSp-ul6a slow', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
setCurrentFigureSize(937, 420)

%% SSp-ul6a fast contra
idxArea = 'SSp-ul6a';
psthPRE = [];
psthPOST = [];
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPRE = [psthPRE m];
            [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthPOST = [psthPOST m];
            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_PRE)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_PRE = [auROC_PRE; auROC];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_POST)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
            end
            [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
            auROC_POST = [auROC_POST; auROC];
        end
    end
end

zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
xlabel('time (s)')
ylabel('zscored firing rate')
[h p] = ttest(auROC_PRE, auROC_POST)
nexttile
plot(auROC_PRE, auROC_POST, 'o')
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'SSp-ul6a fast', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
setCurrentFigureSize(937, 420)


%% striatum tactile fast contra
idxArea = 'CP';
psthPRE = [];
psthPOST = [];
auROC_PRE = [];
auROC_POST = [];
auROC_delta = [];
significant_delta = [];
tot_units = 0;

for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area

    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        ssON = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.phase == 2) & (timestamps.speeds == 2));
        ssOFF = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.phase == 2) & (timestamps.speeds == 2));
        
        for idxUnit = 1:size(st,1) %units in the striatum for each mouse
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, ssOFF, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            countsBsl = [];
            countsRsp = [];

            for idxTrial = 1:length(ssOFF)
                countsBsl(idxTrial) = CountInIntervals(spikes,[ssON(idxTrial)-1 ssON(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[ssOFF(idxTrial) ssOFF(idxTrial)+0.5]);
            end
            [p h] = signrank(countsBsl, countsRsp);
            if h == 1 & countsRsp > countsBsl
                [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPRE = [psthPRE m];
                [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPOST = [psthPOST m];

                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_PRE)
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
                    deltaPRE(idxTrial) = countsRsp(idxTrial) - countsBsl(idxTrial);
                end
                [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
                auROC_PRE = [auROC_PRE; auROC];


                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_POST)
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
                    deltaPOST(idxTrial) = countsRsp(idxTrial) - countsBsl(idxTrial);
                end
                [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
                auROC_POST = [auROC_POST; auROC];

                [auROC, significant] = findAuROC(deltaPRE, deltaPOST, 1, 100);
                auROC_delta = [auROC_delta auROC];
                significant_delta = [significant_delta significant];
                tot_units = tot_units + 1;
            end
        end
    end
end



zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
xlabel('time (s)')
ylabel('zscored firing rate')
[h p] = ttest(auROC_PRE, auROC_POST)
nexttile
plot(auROC_PRE, auROC_POST, 'o')
xlim([0.3 1]);
ylim([0.3 1]);
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'CP fast', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
%setCurrentFigureSize(937, 420)

p_leanrn_units = (sum(significant_delta) ./ tot_units) * 100

%% SSp-ul6a tactile fast contra
idxArea = 'SSp-ul6a';
psthPRE = [];
psthPOST = [];
auROC_PRE = [];
auROC_POST = [];
auROC_delta = [];
significant_delta = [];
tot_units = 0;
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        ssON = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.phase == 2) & (timestamps.speeds == 2));
        ssOFF = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.phase == 2) & (timestamps.speeds == 2));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, ssOFF, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(ssOFF)
                countsBsl(idxTrial) = CountInIntervals(spikes,[ssON(idxTrial)-1 ssON(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[ssOFF(idxTrial) ssOFF(idxTrial)+0.5]);
            end
            [p h] = signrank(countsBsl, countsRsp);
            if h == 1 & countsRsp > countsBsl
                [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPRE = [psthPRE m];
                [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPOST = [psthPOST m];

                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_PRE)
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
                    deltaPRE(idxTrial) = countsRsp(idxTrial) - countsBsl(idxTrial);
                end
                [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
                auROC_PRE = [auROC_PRE; auROC];


                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_POST)
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
                    deltaPOST(idxTrial) = countsRsp(idxTrial) - countsBsl(idxTrial);
                end
                [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
                auROC_POST = [auROC_POST; auROC];

                [auROC, significant] = findAuROC(deltaPRE, deltaPOST, 1, 100);
                auROC_delta = [auROC_delta auROC];
                significant_delta = [significant_delta significant];
                tot_units = tot_units + 1;
            end
        end
    end
end



zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
xlabel('time (s)')
ylabel('zscored firing rate')
[h p] = ttest(auROC_PRE, auROC_POST)
nexttile
plot(auROC_PRE, auROC_POST, 'o')
xlim([0.3 1]);
ylim([0.3 1]);
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'SSp-ul6a fast', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
setCurrentFigureSize(937, 420)

p_leanrn_units = (sum(significant_delta) ./ tot_units) * 100

%% striatum tactile fast contra
idxArea = 'CP';
psthPRE = [];
psthPOST = [];
auROC_PRE = [];
auROC_POST = [];
auROC_delta = [];
significant_delta = [];
tot_units = 0;
for idxMouse  = 1:size(spikeTimes,2)
    timestamps = eventTimes(idxMouse).mouse; %event table of an individual mouse
    mouseT = spikeTimes(idxMouse).mouse; %spikes table of an individual mouse
    st = mouseT{matches(mouseT.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area
    if ~isempty(st)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 1) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.phase == 3) & (timestamps.speeds == 2) & (timestamps.isDarkTrial == 0));
        ssON = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.phase == 2) & (timestamps.speeds == 2));
        ssOFF = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.phase == 2) & (timestamps.speeds == 2));
        for idxUnit = 1:size(st,1)
            spikes = st{idxUnit};
            [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_PRE)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
            end
            [p h] = signrank(countsBsl, countsRsp);
            hPRE = h;

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(visON_POST)
                countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
            end
            [p h] = signrank(countsBsl, countsRsp);
            hPOST= h;

            if hPRE == 1 | hPOST == 1
                [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);     % compute spike raster data
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPRE = [psthPRE m];
                [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPOST = [psthPOST m];

                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_PRE)
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-1 visON_PRE(idxTrial)-0.5]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+0.5]);
                    deltaPRE(idxTrial) = countsRsp(idxTrial) - countsBsl(idxTrial);
                end
                [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
                auROC_PRE = [auROC_PRE; auROC];


                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_POST)
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-1 visON_POST(idxTrial)-0.5]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+0.5]);
                    deltaPOST(idxTrial) = countsRsp(idxTrial) - countsBsl(idxTrial);
                end
                [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
                auROC_POST = [auROC_POST; auROC];

                [auROC, significant] = findAuROC(deltaPRE, deltaPOST, 1, 100);
                auROC_delta = [auROC_delta auROC];
                significant_delta = [significant_delta significant];
                tot_units = tot_units + 1;
            end
        end
    end
end



zscored_psthPRE = zscore(psthPRE, 0, 1);
zscored_psthPOST = zscore(psthPOST, 0, 1);

figure;
L = tiledlayout(1,2);
nexttile
plot(t, mean(zscored_psthPRE'));
hold on
plot(t, mean(zscored_psthPOST'));
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
legend('PRE', 'POST')
xlabel('time (s)')
ylabel('zscored firing rate')
[h p] = ttest(auROC_PRE, auROC_POST)
nexttile
plot(auROC_PRE, auROC_POST, 'o')
xlim([0.3 1]);
ylim([0.3 1]);
stringa = sprintf('p = %0.2f', p);
title(stringa)
xlabel('auROC PRE')
ylabel('auROC POST')
axis square
title(L,'CP fast', 'FontSize', 16, 'Color', 'Red')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
%setCurrentFigureSize(937, 420)

p_leanrn_units = (sum(significant_delta) ./ tot_units) * 100