% This script is the code used to NCT analyse the HCP MEG subjects 
% Compiled by Aradia Wilms, based on scripts from Alina Podschun &
% Sebastian Markett, adapted by Chiara Werner

% starting point: cifti-resampled files for all (83) subjects 

% requirements: this needs the following scripts:
% - firstlevel_ROI_AW.m
% - region_reordering_AW.m 
% - nct_analysis_task_AW.m
% - code by Urs Braun: network_control_and_dopamine-main (Note: the path to
% this script needs to be specified in nct_analysis_task_AW.m)


% files required:
% - functional MEG data in 4k resolution (pre parcellation)
% - an individual dlabel.nii file in 4k resolution 
% - the CATO order aparc.mat



%% 00: setting up

restoredefaultpath;
rehash toolboxcache;


% I chose not to use cd or addpath here for coding efficiency
% so, here we are setting up the general folder stems

% data
data_path = '/home/wernechi/NCT/data/';

% scripts
% we want to add this path for simplicity when calling the functions later
% and it is only one path compared to all the numerous subject paths
addpath '/home/wernechi/NCT/scripts'
addpath '/home/wernechi/cifti-matlab-master' 


%results
results_path = '/home/wernechi/NCT/results/';

% CATO files are not stored in my personal account
CATO_path = '/slow/projects/HCP_1200/01_complete_batch/';


% filenames have the format: 
% '100307_MEG_Wrkmem_srcavgdics_[LM-TIM-0B]_[FB-theta].power.dtseries.nii';


% subject list
sub_string = fileread('subjects.txt');
subj_ids = split(sub_string);
subj_ids = subj_ids(~cellfun('isempty', subj_ids));


% create the inserts for the filenames for later
conditions = {'0B'; '2B'};
freq_bands = {'theta'; 'gammamid'; 'gammalow'; 'gammahigh'; 'delta'; 'betalow'; 'betahigh'; 'alpha'};


%% 01: extract first-level ROIs & region_reordering


for subj = 1:numel(subj_ids)
    
    % open the dlabel file
    atlas_filename = fullfile(data_path, subj_ids{subj}, 'resampled',[subj_ids{subj} '.4k.dlabel.nii']);
    disp(['Aktuelles Subjekt: ', subj_ids{subj}]);
    disp(['Erstelle Pfad: ', atlas_filename]);
    atlas = ciftiopen(atlas_filename);

    % directory in data where this step should be saved
    path = [subj_ids{subj} '/ROIs_and_reordered'];


    for cond = 1:numel(conditions)
        for freq = 1:numel(freq_bands)

            % open the current functional file
            func_file = [subj_ids{subj} '_MEG_Wrkmem_srcavgdics_[LM-TIM-' conditions{cond} ']_[FB-' freq_bands{freq} '].power.dtseries.nii'];
            func_file_path = fullfile(data_path, subj_ids{subj}, func_file);
            cii = ciftiopen(func_file_path);


            % STEP 1: extracting first level ROIs

            % carry out extraction
            [mean_firstlevel]=firstlevel_ROI(atlas, cii);

            % filename
            name = ['mean_firstlevel_ROI_' subj_ids{subj} '_' conditions{cond} '_' freq_bands{freq} '.mat'];
            filename = fullfile(data_path, path, name);
	    folder = fullfile(data_path, path);  % Pfad zum Ordner
	    if ~exist(folder, 'dir')             % Prüfen, ob Ordner existiert
		    mkdir(folder)                % Ordner erstellen, falls nicht vorhanden
	    end



            % due to partial running of this script some participants will
            % have been analysed already: this ensures they are skipped and
            % no errors occur
            if isfile(filename)
            
                fprintf('File %s already exists for %s. Skipping...\n', filename, subj_ids{subj});
                continue; % Skip to the next person in the loop
            else
                try
                    save(filename, 'mean_firstlevel');
                catch ME
                    disp(['Failed to save file', ME.message])
                end
            end
            
            
            % STEP 2: reorder to match CATO parcellation

            CATO_full = fullfile([CATO_path subj_ids{subj} '/DWI_processed_v311/' subj_ids{subj} '_connectivity_gqi_dti_aparc.mat']);
            
            % check whether this person has a CATO file - not all subjects
            % in the wrkmem group have CATO files
            if isfile(CATO_full)
                % load CATO file
                load (CATO_full);
                goal = regionDescriptions;


                % load the first level ROIs
                load (filename);
                eval(['inp = mean_firstlevel;'])


                % do the actual region reordering            
                [data_reord]=region_reordering(goal, atlas, inp);

            
                reord_file = [subj_ids{subj} '_' conditions{cond} '_' freq_bands{freq} '_func_reordered.mat'];
                reord_filename = fullfile(data_path, path, reord_file);


                % skip if file already exists
                if isfile(reord_filename)
                    fprintf('File %s already exists for %s. Skipping...\n', reord_filename, subj_ids{subj});
                    continue; % Skip to the next person in the loop
                else
                    try
                        save(reord_filename,'data_reord');
                    catch ME
                        disp(['HCP ' subj_ids{subj} ' failed to save file', ME.message])
                    end
                end

            else
                continue

            end
            



            % clear out all variables for security
            clear func_file
            clear name
            clear filename
            clear cii
            clear mean_firstlevel
            clear goal
            clear inp
            clear data_reord
            clear reord_file
            clear reord_filename
            clear CATO_full


        end
    end

    clear atlas
    clear atlas_filename
    clear path

    fprintf('Subject with ID %s has been completed.\n', subj_ids{subj});

end


clear subj
clear cond
clear freq

%% 02: actual NCT analysis

% 02.01: prepare input matrices for and then calculate control energy

struct = zeros(68,68,length(subj_ids));

% here I am MANUALLY removing the three persons without a CATO file - this
% makes running the loop easier 
remove_ids = {'662551','102816','169040'};
subj_ids(ismember(subj_ids, remove_ids)) = [];

for subj=1:numel(subj_ids)

    % load CATO order mat
    CATO_load = fullfile([CATO_path subj_ids{subj} '/DWI_processed_v311/' subj_ids{subj} '_connectivity_gqi_dti_aparc.mat']);
    load (CATO_load);
    fprintf('Loaded CATO file for %s \n', subj_ids{subj});

    % this moves FA (in 3) for every participant to the struc matrix
    % CAREFUL! connectivity still contains the subcortical areas
    % these are 14 regions: so we start at the 15th row and column
    struct(:,:,subj) = connectivity(15:end,15:end,3);

    %stabilize like this: A_star=A./(eigs(A,1)+1)-eye(size(A,1));
    A_star(:,:,subj) = struct(:,:,subj)./(eigs(struct(:,:,subj),1)+1)-eye(size(struct(:,:,subj),1));
    A = A_star;

end

clear subj

% 02.02: x0 and xf
% in the end we want to compare the results across different frequency
% bands
% so, this section computes the nct analysis separately for each band

T = 1;
rho = 1;


for freq = 1:numel(freq_bands)

    x0 = zeros(68,length(subj_ids));
    xf = zeros(68,length(subj_ids));

    for subj=1:numel(subj_ids)
    
        % 02.02.a: x0 (0-back)

        reord_file01 = [subj_ids{subj} '/ROIs_and_reordered/' subj_ids{subj} '_' conditions{1} '_' freq_bands{freq} '_func_reordered.mat'];
        reord_filename01 = fullfile(data_path, reord_file01);
        fprintf('Loaded zero back for %s \n', subj_ids{subj})
        

        load(reord_filename01);

        x0(:,subj) = data_reord(:,:);

        
        % 02.02.b: xf (2-back)
        reord_file02 = [subj_ids{subj} '/ROIs_and_reordered/' subj_ids{subj} '_' conditions{2} '_' freq_bands{freq} '_func_reordered.mat'];
        reord_filename02 = fullfile(data_path, reord_file02);
        load(reord_filename02);
        fprintf('Loaded two back for %s \n', subj_ids{subj})
        
        xf(:,subj) = data_reord(:,:);


    end


    % 02.03: access nct_analysis_task_AW.m

    [nct_regions nct_global]=nct_analysis_task(A,T,rho,x0,xf);

    filename = ['NCT_0B_2B_' freq_bands{freq} '.mat'];
    nct_file = fullfile([results_path '00_NCT-raw/' filename]);
    
    try
        save(nct_file,['nct_regions'], ['nct_global']);
    catch ME
        disp(['failed to save NCT', ME.message])
    end

    clear x0
    clear xf
    clear nct_regions
    clear nct_global
    clear filename
    clear nct_file
end

clear freq
clear subj



%% 03: post nct transformations

t = 1001;


% here I am MANUALLY removing the three persons without a CATO file - this
% makes running the loop easier 
remove_ids = {'662551','102816','169040'};
subj_ids(ismember(subj_ids, remove_ids)) = [];

% preallocate for regional measures
wm_reg_T1 = zeros(length(subj_ids),68,4);


for freq = 1:numel(freq_bands)

    nct_file = fullfile([results_path '00_NCT-raw/' ['NCT_0B_2B_' freq_bands{freq} '.mat']]);
    load(nct_file);


    % 03.01: global
    for i=1:numel(subj_ids) % loop over subjects

        wm_glob_T1(i,1,1:2)=1./(sum(nct_global{i}.U(:,3:4))/t); 
        % stability; careful: in original code, 3 and 4 is stability but we want it in 1 and 2. 

        wm_glob_T1(i,1,3:4)=sum(nct_global{i}.U(:,1:2))/t; 
        % control energy; check position as well, respectively
    end


    % 03.02: regional
   
    for s = 1:length(subj_ids) % loop over subjects
        for x = 1:68 % loop over roi
            wm_reg_T1(s,x,1) = 1./(sum(nct_regions{s}.U1(:,x))/t); % stability state A
            wm_reg_T1(s,x,2) = 1./(sum(nct_regions{s}.U2(:,x))/t); % stability state B
            wm_reg_T1(s,x,3) = sum(nct_regions{s}.U3(:,x))/t; % energy A->B
            wm_reg_T1(s,x,4) = sum(nct_regions{s}.U4(:,x))/t; % energy B->A

        end
    end

    % 03.03 saving
    filename1 = fullfile([results_path '01_NCT-final/' freq_bands{freq} '_wm_0B_2B_reg_T1.mat']);
    save(filename1, 'wm_reg_T1');


    filename2 = fullfile([results_path '01_NCT-final/' freq_bands{freq} '_wm_0B_2B_glob_T1.mat']);
    save(filename2, 'wm_glob_T1');


    clear filename1
    clear filename2
    clear wm_glob_T1
    clear wm_reg_T1


end
