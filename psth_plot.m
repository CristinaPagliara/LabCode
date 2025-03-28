%% PSTH without GUI

%% 1 import and set
expRootFolder = 'E:\PhD\E04_PERIPE_PP-CP\Data\explorerGUI_dataset';
load(fullfile(expRootFolder, 'pooledTable.mat'), 'T')
load(fullfile(expRootFolder, 'pooledDataSet.mat'), 'ts')

%% 1.1 extract selected units and trials 
units_conditions=  T.pANOVA_baseVStest>0.05;%T.name == "Caudoputamen" &
T_units=find(units_conditions)';

for n=1:sum(units_conditions)
    neuron=T_units(n);
    cellID=T.uID(T_units(n));
    ms=T.mouseNum(neuron);
    for t=1:lenght(ts(ms).T)
        trials_conditions= ts(ms).T.isLit==1 & ts(ms).T.isFar == 1; %COME SELEZIONO IL BLOCCO IN MANIERA UNIVOCA PER TUTTI I TOPI????
    end
end

%% 1.1 extract selected trials: change the neuron variable and the trial conditions
neuron = 7;

mn = T.mouseNum(neuron);
cuse = T.uID(neuron);

trialsuse=ts(mn).T.isLit==1 & ts(mn).T.isFar == 1; %the number of trials that satisfied that conditions
trialnumbers = find(trialsuse)'; %all trials that satisfied the trialuse conditions in mouse T_ordinal   

%% 1.2 compute psth
psthTiming.preWindow = 2;
psthTiming.postWindow = 3;
psthTiming.binSize = 0.05;
psthTimimg.baselineWindow = [-2, 0];  
edges = -psthTiming.preWindow:psthTiming.binSize:psthTiming.postWindow;
smoothingw = 1;

thisbaselinespikerate = nan(sum(trialsuse),1);

for idx = 1:sum(trialsuse)
    ntr = trialnumbers(idx);
    spiketimes = ts(mn).T.raster{ntr};
    cellIDs = ts(mn).T.cellIDs{ntr};
    spiketimes = spiketimes(cellIDs==cuse);
    thisbaselinespikerate(ntr) = sum( spiketimes>psthTimimg.baselineWindow(1) & spiketimes<=psthTimimg.baselineWindow(2) ) / diff(psthTimimg.baselineWindow);
    % baseline is always calculated before the first event
    psth_trial = histcounts(spiketimes, edges)./psthTiming.binSize - thisbaselinespikerate(ntr);
    psth(idx,:) = gen_fx_gsmooth(psth_trial, smoothingw);
end

%% 1.3 display psth
psth_time = edges(1:end-1) + psthTiming.binSize/2;
baseColor = [0.5, 0.2, 0.45];

fapp = figure('Visible', 'on');

% main psth plot with error shade:
stdshade_modified(psth,0.35,baseColor,psth_time,[],[], 2.5, 2);

% Formatting (removes axes and labels)
axis off; % Completely remove axes
xline(0, 'Color', [0.3804    0.1059    0.2902], 'LineStyle', '--', 'LineWidth', 1.5);




























