clear
clc

expRootFolder = 'E:\PhD\E04_PERIPE_PP-CP\Data\ephys analyses';
load(fullfile(expRootFolder, 'pooledTable.mat'), 'T')
load(fullfile(expRootFolder, 'pooledDataSet.mat'), 'ts')

for idxMouse = 1:size(ts,2)
    mouse(idxMouse).eventTimes = ts(idxMouse).T(:, ["sideIdx", "endloc", "speeds", "phase", "isDarkTrial", "Ton", "Toff", "Ton_rec", "Toff_rec"]);
    mouse(idxMouse).spikeTimes = T(T.mouseNum == idxMouse, ["kID", "acronym", "templWaveform", "spikeTimes"]);
end
clear T ts

%% 1 striatum tactile fast contra - among the tactile units how many have a visual potentiation? 
idxArea = 'CP';

psthCONTACT=[];
psthPRE=[];
psthPOST=[];

auROC_PRE=[];
auROC_POST=[];
auROC_delta = [];
significant_delta=[];

tactile_units=[];

tot_units=0;

for idxMouse  = 1:size(mouse,2)
    timestamps = mouse(idxMouse).eventTimes; %event table of an individual mouse
    mouse_units = mouse(idxMouse).spikeTimes; %spikes table of an individual mouse
    area_spikes = mouse_units{matches(mouse_units.acronym, idxArea), 'spikeTimes'}; %spiketimes from a specific area

    if ~isempty(area_spikes)
        visON_PRE = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.speeds == 2) & (timestamps.phase == 1) & (timestamps.isDarkTrial == 0));
        visOFF_PRE = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1)  & (timestamps.speeds == 2) & (timestamps.phase == 1) & (timestamps.isDarkTrial == 0));
        visON_POST = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.speeds == 2) & (timestamps.phase == 3) & (timestamps.isDarkTrial == 0));
        visOFF_POST = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 1) & (timestamps.speeds == 2) & (timestamps.phase == 3) & (timestamps.isDarkTrial == 0));
        ssON = timestamps.Ton((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.speeds == 2) & (timestamps.phase == 2));
        ssOFF = timestamps.Toff((timestamps.sideIdx == 1) & (timestamps.endloc == 3) & (timestamps.speeds == 2) & (timestamps.phase == 2));
        
        sig_contact=[]; 

        for unit = 1:size(area_spikes,1) %units in the IdxArea for each mouse
            spikes = area_spikes{unit};

            [raster,indices] = Sync(spikes, ssOFF, 'durations', [-1 1]);
            [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
            psthCONTACT = [psthCONTACT m];

            countsBsl = [];
            countsRsp = [];
            for idxTrial = 1:length(ssOFF)
                % Trial_number= find(ismember(timestamps.Toff, ssOFF(idxTrial)));
                window = min(timestamps.Ton_rec-timestamps.Toff);
                countsBsl(idxTrial) = CountInIntervals(spikes,[ssON(idxTrial)-window ssON(idxTrial)]);
                countsRsp(idxTrial) = CountInIntervals(spikes,[ssOFF(idxTrial) ssOFF(idxTrial)+window]);
            end
            [p h] = signrank(countsBsl, countsRsp);
            sig_contact= [sig_contact h];

            if h == 1 
                [raster,indices] = Sync(spikes, visON_PRE, 'durations', [-1 1]);    
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPRE = [psthPRE m];

                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_PRE)
                    window = min(timestamps.Toff-timestamps.Ton);
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial)-window visON_PRE(idxTrial)]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_PRE(idxTrial) visON_PRE(idxTrial)+window]);
                    deltaPRE(idxTrial) = countsRsp(idxTrial) - countsBsl(idxTrial);
                end
                [auROC, significant] = findAuROC(countsBsl, countsRsp, 0, 100);
                auROC_PRE = [auROC_PRE; auROC];

                [raster,indices] = Sync(spikes, visON_POST, 'durations', [-1 1]);     % compute spike raster data
                [m,t] = SyncHist(raster, indices,'mode','mean','durations', [-1 1], 'smooth', 1);
                psthPOST = [psthPOST m];

                countsBsl = [];
                countsRsp = [];
                for idxTrial = 1:length(visON_POST)
                    window = min(timestamps.Toff-timestamps.Ton);
                    countsBsl(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial)-window visON_POST(idxTrial)]);
                    countsRsp(idxTrial) = CountInIntervals(spikes,[visON_POST(idxTrial) visON_POST(idxTrial)+window]);
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
        tactile_units(1,idxMouse)=sum(sig_contact);
    end
end

p_leanrn_units = (sum(significant_delta) ./ tot_units) * 100

%
figure;

subplot(1,3,1);
zscored_psthCONTACT = zscore(psthCONTACT, 0, 1);
stdshade_modified(mean(zscored_psthCONTACT'),0.35,[0.5, 0.2, 0.45],t,[],[], 2.5, 2);
title('contact');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
xlabel('time (s)')
ylabel('zscored firing rate')

subplot(1,3,2);
zscored_psthPRE = zscore(psthPRE, 0, 1);
stdshade_modified(mean(zscored_psthPRE'),0.35,[0.4, 0.7, 0.7],t,[],[], 2.5, 2);
title('pre');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
xlabel('time (s)')
ylabel('zscored firing rate')

subplot(1,3,3);
zscored_psthPOST = zscore(psthPOST, 0, 1);
stdshade_modified(mean(zscored_psthPOST'),0.35,[0.075, 0.400, 0.510],t,[],[], 2.5, 2);
title('post');

set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
xlabel('time (s)')
ylabel('zscored firing rate')

%% 2 striatum visual fast contra - among the visual units how many have a visual potentiation? 