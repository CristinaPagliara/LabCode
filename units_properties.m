clear
clc

% 1 import and set
expRootFolder = 'E:\PhD\E04_PERIPE_PP-CP\Data\ephys analyses';
load(fullfile(expRootFolder, 'pooledTable.mat'), 'T')
load(fullfile(expRootFolder, 'pooledDataSet.mat'), 'ts')

%% 1.1 extract selected units and trials 
T_units=find(T.name == "Caudoputamen")';
numUnits = length(T_units);

%fix this 
contact_firing=nan(58,numUnits);
hit_baseline=nan(58,numUnits);
app_firing=nan(152,numUnits);
far_baseline=nan(152,numUnits);
rec_firing=nan(152,numUnits);

meanFiring_hit = nan(numUnits, 2);
meanFiring_far = nan(numUnits, 3);

p_values_hit = nan(1, numUnits);
p_values_app = nan(1, numUnits);
p_values_rec = nan(1, numUnits);

for n=1:numUnits %loop on caudoputamen units
    neuron=T_units(n);
    ms=T.mouseNum(neuron);
    cellID=T.uID(T_units(n));

    hit_trials= find(ts(ms).T.endloc == 3 & ts(ms).T.phase == 2 & ts(ms).T.sideIdx == 1 & ts(ms).T.speeds == 2)';
    % near_trials= find(ts(ms).T.isNear == 1)';
    % contact_trials= find(ts(ms).T.isNear == 1 | ts(ms).T.isHit == 1)';
    far_trials= find(ts(ms).T.endloc == 1 & ts(ms).T.isLit == 1 & ts(ms).T.isDarkTrial == 0 & ts(ms).T.sideIdx == 1 & ts(ms).T.speeds == 2)';

    for hit_t=1:length(hit_trials) %loop on hit_trials
        trial_hit=hit_trials(hit_t);
        %hit_baseline(hit_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Toff_rec(trial_hit-1) & T.spikeTimes{neuron,1}<ts(ms).T.Ton(trial_hit))/(ts(ms).T.Ton(trial_hit)-ts(ms).T.Toff_rec(trial_hit-1));
        %contact_firing(hit_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Toff(trial_hit) & T.spikeTimes{neuron,1}<ts(ms).T.Ton_rec(trial_hit))/(ts(ms).T.Ton_rec(trial_hit)-ts(ms).T.Toff(trial_hit));
        hit_baseline(hit_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Ton(trial_hit)-1 & T.spikeTimes{neuron,1}<ts(ms).T.Ton(trial_hit)-0.5)/((ts(ms).T.Ton(trial_hit)-0.5) - (ts(ms).T.Ton(trial_hit)-1));
        contact_firing(hit_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Toff(trial_hit) & T.spikeTimes{neuron,1}<ts(ms).T.Toff(trial_hit)+0.5)/((ts(ms).T.Toff(trial_hit)+0.5) -ts(ms).T.Toff(trial_hit));        
    end

    % Compute the average firing rate over all trials for each condition
    meanFiring_hit(n, 1) = nanmean(hit_baseline(:, n));
    meanFiring_hit(n, 2) = nanmean(contact_firing(:, n));

    % Run the Wilcoxon signed-rank test for this neuron over all trials
    p = signrank(hit_baseline(:, n), contact_firing(:, n),'method', 'exact');
    p_values_hit(n) = p;

    if p_values_hit(n) < 0.05 % amonge the tactile units
        for far_t=1:length(far_trials) %loop on the
            trial_far=far_trials(far_t);
            if trial_far>1
                far_baseline(far_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Toff_rec(trial_far-1) & T.spikeTimes{neuron,1}<ts(ms).T.Ton(trial_far))/(ts(ms).T.Ton(trial_far)-ts(ms).T.Toff_rec(trial_far-1));
                app_firing(far_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Ton(trial_far) & T.spikeTimes{neuron,1}<ts(ms).T.Toff(trial_far))/(ts(ms).T.Toff(trial_far)-ts(ms).T.Ton(trial_far));
                rec_firing(far_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Ton_rec(trial_far) & T.spikeTimes{neuron,1}<ts(ms).T.Toff_rec(trial_far))/(ts(ms).T.Toff_rec(trial_far)-ts(ms).T.Ton_rec(trial_far));
            end
            if trial_far==1
                far_baseline(far_t,n)=sum(T.spikeTimes{neuron,1}<ts(ms).T.Ton(trial_far))/(ts(ms).T.Ton(trial_far)-T.spikeTimes{neuron,1}(1,1));
                app_firing(far_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Ton(trial_far) & T.spikeTimes{neuron,1}<ts(ms).T.Toff(trial_far))/(ts(ms).T.Toff(trial_far)-ts(ms).T.Ton(trial_far));
                rec_firing(far_t,n)=sum(T.spikeTimes{neuron,1}>ts(ms).T.Ton_rec(trial_far) & T.spikeTimes{neuron,1}<ts(ms).T.Toff_rec(trial_far))/(ts(ms).T.Toff_rec(trial_far)-ts(ms).T.Ton_rec(trial_far));
            end
        end
    else
        far_baseline(1:length(far_trials), n) = NaN;
        app_firing(1:length(far_trials), n) = NaN;
        rec_firing(1:length(far_trials), n) = NaN;
    end
    
    % Check if far_baseline is empty before running signrank
    if ~all(isnan(far_baseline(:, n)))
        pa= signrank(far_baseline(:, n), app_firing(:, n), 'method', 'exact');
        p_values_app(n) = pa;
        pr = signrank(far_baseline(:, n), rec_firing(:, n), 'method', 'exact');
        p_values_rec(n) = pr;
        meanFiring_far(n, 1) = nanmean(far_baseline(:, n));
        meanFiring_far(n, 2) = nanmean(app_firing(:, n));
        meanFiring_far(n, 3) = nanmean(rec_firing(:, n));
    else
        p_values_app(1,n) = NaN;
        p_values_rec(1,n) = NaN;
        meanFiring_far(n,1:3)= NaN;
    end
end

%%
% Identify indices for significant and non-significant neurons
sigIdx_hit = p_values_hit < 0.05;
nonSigIdx_hit = ~sigIdx_hit;

sigIdx_app = p_values_app < 0.05;
nonSigIdx_app = ~sigIdx_app;

sigIdx_rec = p_values_rec < 0.05;
nonSigIdx_rec = ~sigIdx_rec;

sigIdx_far = sigIdx_rec+ sigIdx_app;

% figure;
% scatter(meanFiring_hit(nonSigIdx_hit, 1), meanFiring_hit(nonSigIdx_hit, 2), 'filled', 'MarkerFaceColor', [153,204,255]/255);
% hold on;
% scatter(meanFiring_hit(sigIdx_hit, 1), meanFiring_hit(sigIdx_hit, 2), 'filled', 'MarkerFaceColor', [204,0,102]/255);
% 
% xlabel('Baseline Firing Rate (Hz)');
% ylabel('Contact Firing Rate (Hz)');
% title('Baseline vs. Contact Firing Rates for Each Neuron');
% 
% % Plot a 45° reference line for equality
% lims = [min([meanFiring_hit(:); 0]) max(meanFiring_hit(:))];
% plot(lims, lims, '--k');
% 
% %legend('p \geq 0.05', 'p < 0.05', 'Equal Firing Rate', 'Location', 'best');
% 
% % Remove axis ticks, box, and set transparent background
% set(gca, 'XTick', [], 'YTick', [], 'Box', 'off', 'Color', 'none');  
% 
% hold off;
% 
% %%
% % Identify indices for significant and non-significant neurons
% p_values_far = p_values_far(~isnan(p_values_far));
% valid_rows = ~all(isnan(meanFiring_far), 2); 
% meanFiring_far = meanFiring_far(valid_rows, 1:2);
% 
% sigIdx_far = p_values_far < 0.05;
% nonSigIdx_far = ~sigIdx_far;
% 
% figure;
% scatter(meanFiring_far(nonSigIdx_far, 1), meanFiring_far(nonSigIdx_far, 2), 'filled', 'MarkerFaceColor', [153,204,255]/255);
% hold on;
% scatter(meanFiring_far(sigIdx_far, 1), meanFiring_far(sigIdx_far, 2), 'filled', 'MarkerFaceColor', [204,0,102]/255);
% 
% xlabel('Baseline Firing Rate (Hz)');
% ylabel('Contact Firing Rate (Hz)');
% title('Baseline vs. Contact Firing Rates for Each Neuron');
% 
% % Plot a 45° reference line for equality
% lims = [min([meanFiring_far(:); 0]) max(meanFiring_far(:))];
% plot(lims, lims, '--k');
% 
% %legend('p \geq 0.05', 'p < 0.05', 'Equal Firing Rate', 'Location', 'best');
% 
% % Remove axis ticks, box, and set transparent background
% set(gca, 'XTick', [], 'YTick', [], 'Box', 'off', 'Color', 'none');  
% 
% hold off;

%%
sig_contact = find(sigIdx_hit);
sig_contact_app = find(sigIdx_app);
sig_contact_rec = find(sigIdx_rec);

not_sig= find(sigIdx_far(:) == 0);
sig_contact_far = find(sigIdx_far(:) == 2);
sig_contact_far_one = find(sigIdx_far(:) == 1);

numUnits = length(p_values_hit);

% Count the number of units per category
counts = [numUnits, length(sig_contact), (length(sig_contact_app)+length(sig_contact_rec)-length(sig_contact_far)), (length(sig_contact)-(length(sig_contact_app)+length(sig_contact_rec)-length(sig_contact_far))), length(not_sig)];

