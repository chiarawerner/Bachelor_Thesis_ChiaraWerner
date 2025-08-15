% script does anovas for control energy for a task
% based on scripts from Aradia Wilms

%% general stuff; paths etc

addpath '/home/wernechi/NCT/scripts'

%This is the list of subject directory names:
sub_string = fileread('subjects.txt');
subj_ids = split(sub_string);
subj_ids = subj_ids(~cellfun('isempty', subj_ids));

% remove three subj that didn't have a CATO
remove_ids = {'662551','102816','169040'};
subj_ids(ismember(subj_ids, remove_ids)) = [];


% create the inserts for the filenames for later
freq_bands = {'theta'; 'gammamid'; 'gammalow'; 'gammahigh'; 'delta'; 'betalow'; 'betahigh'; 'alpha'};

% NCT measures to test across
measures = {'glob'; 'reg'};



%% prepare mean difference


% compute mean difference for each individual
for s=1:numel(subj_ids)
    for f=1:numel(freq_bands)
        
        load (['/home/wernechi/NCT/data/' subj_ids{s} '/ROIs_and_reordered/' subj_ids{s} '_2B_' freq_bands{f} '_func_reordered.mat']);
        two_back = data_reord;


        load (['/home/wernechi/NCT/data/' subj_ids{s} '/ROIs_and_reordered/' subj_ids{s} '_0B_' freq_bands{f} '_func_reordered.mat']);
        zero_back = data_reord;

        difference_score = zeros(68,1);

        for roi=1:numel(two_back)
            difference_score=two_back - zero_back;
        end
        
        filename=['/home/wernechi/NCT/data/' subj_ids{s} '/ROIs_and_reordered/' subj_ids{s} '_' freq_bands{f} '_mean_dif_activity.mat'];
        save(filename,'difference_score')
    end
end



%% compute average of mean difference
% per person across ROIs for global measures

% Initialize a vector to store the average values for each subject
average_values = zeros(length(subj_ids), 1);

for f = 1:length(freq_bands)
    for s = 1:length(subj_ids)
        % Construct the file path
        file_path = ['/home/wernechi/NCT/data/' subj_ids{s} '/ROIs_and_reordered/' subj_ids{s} '_' freq_bands{f} '_mean_dif_activity.mat'];
        
        % Load the data from the file
        data_struct = load(file_path);
        data = data_struct.difference_score; % Assuming the variable is named 'difference_score'
        
        % Compute the average value across all rows for this subject
        average_value = mean(data, 'all');
        
        % Store the result in the vector
        average_values(s, 1) = average_value;
    end
    filename=['/home/wernechi/NCT/results/02_NCT-tested/00_meandiff/' freq_bands{f} 'mean_dif_activity.mat'];
    save(filename, 'average_values')
end


%% prepare tables and model; do anova


for freq=1:numel(freq_bands)

    load (['/home/wernechi/NCT/results/01_NCT-final/' freq_bands{freq} '_wm_0B_2B_glob_T1.mat']);
    nct_results = wm_glob_T1;

    load(['/home/wernechi/NCT/results/02_NCT-tested/00_meandiff/' freq_bands{freq} 'mean_dif_activity.mat']);
    cv = average_values;

    % STABILITY
    % write observations + cv to table
    tt=table(squeeze(nct_results(:,:,1)),squeeze(nct_results(:,:,2)),cv);

    % within-design
    mm=table([1 2]','VariableNames',{'Measurements'});

    % fit rm model to data in tt, treat first two columns as repeated
    % measauremens, and the variable after ~ as between subject factor
    rm = fitrm(tt,'Var1-Var2~cv','WithinDesign',mm);

    % run anova on rm model
    tbl=ranova(rm);

    % extract p-value for quick inspection
    pval_nct=table2array(tbl(1,5));
        
    % get eta-squared as effect size
    eta_nct=tbl.SumSq(1)/(tbl.SumSq(1)+tbl.SumSq(3));
    eval(['results.stability.tbl = tbl']);
    eval(['results.stability.pval = pval_nct']);
    eval(['results.stability.eta = eta_nct']);

    % ENERGY
    % write observations + cv to table
    tt=table(squeeze(nct_results(:,:,3)),squeeze(nct_results(:,:,4)),cv);
    % within-design
    mm=table([1 2]','VariableNames',{'Measurements'});
    % fit rm model to data in tt, treat first two columns as repeated
    % measauremens, and the variable after ~ as between subject factor
    rm = fitrm(tt,'Var1-Var2~cv','WithinDesign',mm);
    % run anova on rm model
    tbl=ranova(rm);
    % extract p-value for quick inspection
    pval_nct=table2array(tbl(1,5));
    % get eta-squared as effect size
    eta_nct=tbl.SumSq(1)/(tbl.SumSq(1)+tbl.SumSq(3));
    eval (['results.energy.tbl = tbl']);
    eval (['results.energy.pval = pval_nct']);
    eval (['results.energy.eta = eta_nct']);


    savenct=['/home/wernechi/NCT/results/02_NCT-tested/' freq_bands{freq} '_stats_glob.mat'];
    save(savenct,['results']);


end

%% post hoc paired t-test 

% --- PAIRED T-TEST: STABILITY ---
stability_0B = squeeze(nct_results(:,:,1));  % 0-back
stability_2B = squeeze(nct_results(:,:,2));  % 2-back

[~, p_stability, ~, stats_stability] = ttest(stability_0B, stability_2B);

% Ergebnisse speichern
results.stability.ttest.tstat = stats_stability.tstat;
results.stability.ttest.df = stats_stability.df;
results.stability.ttest.pval = p_stability;

% --- PAIRED T-TEST: ENERGY ---
energy_0B = squeeze(nct_results(:,:,3));  % 0-back
energy_2B = squeeze(nct_results(:,:,4));  % 2-back

[~, p_energy, ~, stats_energy] = ttest(energy_0B, energy_2B);

% Ergebnisse speichern
results.energy.ttest.tstat = stats_energy.tstat;
results.energy.ttest.df = stats_energy.df;
results.energy.ttest.pval = p_energy;
