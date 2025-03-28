% 1 import and set
expRootFolder = 'E:\PhD\E04_PERIPE_PP-CP\Data\explorerGUI_dataset';
load(fullfile(expRootFolder, 'pooledTable.mat'), 'T')
load(fullfile(expRootFolder, 'pooledDataSet.mat'), 'ts')

%% 1.1 extract selected units and trials 
units_conditions=  T.pANOVA_baseVStest>0.05 & T.name == "Caudoputamen"; %filter units
T_units=find(units_conditions)';



for n=1:length(T_units) %loop among units
    neuron=T_units(n);
    cellID=T.uID(T_units(n));
    ms=T.mouseNum(neuron);

    hit_trials= find(ts(ms).T.isHit == 1)';
    near_trials= find(ts(ms).T.isNear == 1)';
    contact_trials= find(ts(ms).T.isNear == 1 | ts(ms).T.isHit == 1)';
    
    for t=1:length(hit_trials) %loop among trials
        trial=hit_trials(t);
        trial_ID=ts(ms).T.cellIDs{trial,1}==cellID;
        spiketimes=ts(ms).T.raster{trial,1}(trial_ID);
        % what is the 0 in the raster column? 
        %spike_contact_phase=sum(spiketimes>(ts(ms).T.Toff(trial,1)-ts(ms).T.Ton(trial,1)) & spiketimes<(ts(ms).T.Ton_rec(trial,1)-ts(ms).T.Ton(trial,1)))/;
        %spike_baseline = sum(spiketimes<ts(ms));
    end
end
