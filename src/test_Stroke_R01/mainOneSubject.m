%% Created by MT 02/04/25
% The main pipeline for R01 Stroke Spinal Stim Aim 1 (using tables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Comment this part out when running all subjects at once.
% clc;
% clearvars;
% subject = 'SS01';
% configFilePath = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\src\overground\config.json';
% config = jsondecode(fileread(configFilePath));
% disp(['Loaded configuration from: ' configFilePath]);
% doPlot = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get configuration
intervention_folders = config.INTERVENTION_FOLDERS;
intervention_field_names = config.MAPPED_INTERVENTION_FIELDS;
mapped_interventions = containers.Map(intervention_folders, intervention_field_names);
gaitriteConfig = config.GAITRITE;
delsysConfig = config.DELSYS_EMG;
xsensConfig = config.XSENS;
regexsConfig = config.REGEXS;

% Folder to load the data from.
pathsConfig = config.PATHS;
subjectLoadPath = fullfile(pathsConfig.ROOT_LOAD, subject);
% Path to save the data to.
subjectSaveFolder = fullfile(pathsConfig.ROOT_SAVE, subject);
saveFileName = pathsConfig.SAVE_FILE_NAME;
codeFolderPath = pathsConfig.CODE_FOLDER_PATH; % Folder where the code lives
addpath(genpath(pathsConfig.CODE_FOLDER_PATH));

%% Initialize outcome measure tables
trialTable = table; % Each row is one trial, all data
cycleTable = table; % Each row is one UNMATCHED gait cycle, all data
visitTable = table; % Each row is one whole session
speedPrePostTable = table; % Each row is one combination of SSV/FV & Pre/Post
cycleTableContraRemoved = table; % Each row is one UNMATCHED gait cycle, with the contralateral data removed and column names merged

%% GaitRite Load
subject_gaitrite_folder = fullfile(subjectLoadPath, gaitriteConfig.FOLDER_NAME);
gaitRiteTable = loadGaitRiteAllInterventions(gaitriteConfig, subject_gaitrite_folder, intervention_folders, mapped_interventions, regexsConfig);

%% Load Delsys
subject_delsys_folder = fullfile(subjectLoadPath, delsysConfig.FOLDER_NAME);
delsysTable = loadDelsysAllInterventions(delsysConfig, subject_delsys_folder, intervention_folders, mapped_interventions, regexsConfig);

%% Filter Delsys
delsysTableFiltered = filterDelsys(delsysTable, 'Delsys_Loaded', 'Delsys_Filtered', delsysConfig.FILTER, delsysConfig.SAMPLING_FREQUENCY);
delsysTable = addToTable(delsysTable, delsysTableFiltered);

%% Load XSENS
subject_xsens_folder = fullfile(subjectLoadPath, xsensConfig.FOLDER_NAME);
xsensTable = loadXSENSAllInterventions(xsensConfig, subject_xsens_folder, intervention_folders, mapped_interventions, regexsConfig);

%% Filter XSENS
xsensTableFiltered = filterXSENS(xsensTable, 'XSENS_Loaded', 'XSENS_Filtered', xsensConfig.FILTER, xsensConfig.SAMPLING_FREQUENCY);
xsensTable = addToTable(xsensTable, xsensTableFiltered);

%% Correct the GaitRite trial numbers when there's fewer than 3 trials
gaitRiteTable = correctGRTrialNumbersWhenLessThan3(gaitRiteTable, {xsensTable, delsysTable});

%% Adjust the order of trials to match GaitRite, as needed
[xsensTableReordered, delsysTableReordered] = checkTrialOrderAllInterventions(gaitRiteTable, {xsensTable, delsysTable});
trialTable = addToTable(trialTable, gaitRiteTable);
trialTable = addToTable(trialTable, delsysTableReordered);
trialTable = addToTable(trialTable, xsensTableReordered);
dateTimeColNames = trialTable.Properties.VariableNames(contains(trialTable.Properties.VariableNames, 'DateTimeSaved'));
trialTable = removevars(trialTable, dateTimeColNames);

%% Plot raw and filtered timeseries data
% if doPlot
%     baseSavePathEMG = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\EMG\Raw_Filtered';
%     plotRawAndFilteredData(trialTable, 'Raw and Filtered EMG', baseSavePathEMG, struct('Raw','Delsys_Loaded', 'Filtered', 'Delsys_Filtered'), true);
%     baseSavePathXSENS = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\Joint Angles\Raw_Filtered';
%     plotRawAndFilteredData(trialTable, 'Raw and Filtered Joint Angles', baseSavePathXSENS, struct('Raw','XSENS_Loaded', 'Filtered', 'XSENS_Filtered'), false);
% end

%% Time synchronization
syncedTableDelsys = timeSynchronize(trialTable, delsysConfig.SAMPLING_FREQUENCY, 'seconds', 'Delsys_Frames');
syncedTableXSENS = timeSynchronize(trialTable, xsensConfig.SAMPLING_FREQUENCY, 'seconds', 'XSENS_Frames');
trialTable = addToTable(trialTable, syncedTableDelsys);
trialTable = addToTable(trialTable, syncedTableXSENS);

%% Plot each trial's data individually, along with gait event information.
% if doPlot
%     baseSavePathEMG = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\EMG\Trials_GaitEvents';
%     plotTrialWithGaitEvents(trialTable, 'Filtered EMG and Gait Events', baseSavePathEMG, 'Delsys_Filtered', 'Delsys_Frames');
%     baseSavePathXSENS = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\Joint Angles\Trials_GaitEvents';
%     plotTrialWithGaitEvents(trialTable, 'Filtered Joint Angles and GaitEvents', baseSavePathXSENS, 'XSENS_Filtered', 'XSENS_Frames');
% end

%% Split data by gait cycle without doing any matching between L & R gait cycles
xsensCyclesTable = splitTrialsByGaitCycle_NoMatching(trialTable, 'XSENS_Filtered','XSENS_Frames');
delsysCyclesTable = splitTrialsByGaitCycle_NoMatching(trialTable, 'Delsys_Filtered','Delsys_Frames');
cycleTable = addToTable(cycleTable, xsensCyclesTable);
cycleTable = addToTable(cycleTable, delsysCyclesTable);

%% Distribute GaitRite vectors from the trial table to the gait cycle table.
% e.g. step/stride lengths/widths/durations/etc.
% Also include the start and end of each gait cycle and swing/stance phase
grDistributedTable = distributeGaitRiteDataToSeparateTable(gaitRiteTable);

%% Plot each gait cycle's filtered data, non-time normalized and each gait cycle of one condition plotted on top of each other.
% if doPlot
%     baseSavePathEMG = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\EMG\Filtered_GaitCycles';
%     plotAllTrials(delsysCyclesTable, 'Filtered EMG', baseSavePathEMG, 'Delsys_Filtered');
%     baseSavePathXSENS = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\Joint Angles\Filtered_GaitCycles';
%     plotAllTrials(xsensCyclesTable, 'Filtered Joint Angles', baseSavePathXSENS, 'XSENS_Filtered');
% end

%% Downsample each gait cycle's data to 101 points.
n_points = config.NUM_POINTS;
xsensDownsampledTable = downsampleAllData(cycleTable, 'XSENS_Filtered', 'XSENS_TimeNormalized', n_points);
delsysDownsampledTable = downsampleAllData(cycleTable, 'Delsys_Filtered', 'Delsys_TimeNormalized', n_points);
cycleTable = addToTable(cycleTable, xsensDownsampledTable);
cycleTable = addToTable(cycleTable, delsysDownsampledTable);

%% Identify the max EMG data value across one whole visit (all trials & gait cycles)
maxEMGTable = maxEMGValuePerVisit(cycleTable, 'Delsys_TimeNormalized', 'Max_EMG_Value');
visitTable = addToTable(visitTable, maxEMGTable);

%% Normalize the time-normalized EMG data to the max value across one whole visit (all trials & gait cycles)
normalizedEMGTable = normalizeAllDataToVisitValue(cycleTable, 'Delsys_TimeNormalized', visitTable, 'Max_EMG_Value', 'Delsys_Normalized_TimeNormalized', 2);
cycleTable = addToTable(cycleTable, normalizedEMGTable);

%% Get % gait cycle when peaks occur
emgPeaksTable = getPeaks(cycleTable, 'Delsys_Normalized_TimeNormalized', 'EMG', {'max'}, {'index'});
jointAnglesPeaksTable = getPeaks(cycleTable, 'XSENS_TimeNormalized', 'JointAngles', {'min', 'max'}, {'index'});
cycleTable = addToTable(cycleTable, emgPeaksTable);
cycleTable = addToTable(cycleTable, jointAnglesPeaksTable);

%% Plot each gait cycle's scaled to max EMG data, and each gait cycle of one condition plotted on top of each other.
% if doPlot
%     baseSavePathEMG = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\EMG\ScaledToMax_GaitCycles';
%     plotAllTrials(matchedCycleTable, 'Scaled To Max EMG', baseSavePathEMG, 'Delsys_Normalized_TimeNormalized');    
% end

%% Set up muscle & joint names for analyses
disp('Defining L & R names');
musclesLR = delsysConfig.MUSCLES;
jointsLR = xsensConfig.JOINTS;
musclesL = cell(size(musclesLR));
musclesR = cell(size(musclesLR));
jointsL = cell(size(jointsLR));
jointsR = cell(size(jointsLR));
for i = 1:length(musclesLR)
    musclesL{i} = ['L' musclesLR{i}];
    musclesR{i} = ['R' musclesLR{i}];
end
for i = 1:length(jointsLR)
    jointsL{i} = ['L' jointsLR{i}];
    jointsR{i} = ['R' jointsLR{i}];
end

%% Calculate the number of muscle synergies in each gait cycle of each trial
VAFthresh = config.DELSYS_EMG.VAF_THRESHOLD;
synergiesTableL = calculateSynergiesAll(cycleTable, 'Delsys_Normalized_TimeNormalized', musclesL, VAFthresh, 'L');
synergiesTableR = calculateSynergiesAll(cycleTable, 'Delsys_Normalized_TimeNormalized', musclesR, VAFthresh, 'R');
cycleTable = addToTable(cycleTable, synergiesTableL);
cycleTable = addToTable(cycleTable, synergiesTableR);

%% Scatterplot the number of muscle synergies & the step lengths
% if doPlot
%     baseSavePathEMG = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\EMG\NumSynergies';
%     scatterPlotPerGaitCyclePerIntervention(delsysStruct, 'Num Synergies', baseSavePathEMG, 'NumSynergies');    
% end

%% Calculate area under the curve (AUC)
aucTableXSENS = calculateAUCAll(cycleTable, 'XSENS_TimeNormalized', 'AUC_JointAngles');
aucTableDelsys = calculateAUCAll(cycleTable, 'Delsys_Normalized_TimeNormalized', 'AUC_EMG');
cycleTable = addToTable(cycleTable, aucTableXSENS);
cycleTable = addToTable(cycleTable, aucTableDelsys);

%% Calculate root mean square (RMS)
rmsTableXSENS = calculateRMSAll(cycleTable, 'XSENS_TimeNormalized', 'RMS_JointAngles');
rmsTableDelsys = calculateRMSAll(cycleTable, 'Delsys_Normalized_TimeNormalized', 'RMS_EMG');
cycleTable = addToTable(cycleTable, rmsTableXSENS);
cycleTable = addToTable(cycleTable, rmsTableDelsys);

%% Calculate range of motion (ROM)
romTableXSENS = calculateRangeAll(cycleTable, 'XSENS_TimeNormalized', 'JointAngles');
cycleTable = addToTable(cycleTable, romTableXSENS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below this point, computations require the left and right gait cycles to be matched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Match the L and R gait cycles for symmetry analysis
matchedCycleTable = matchCycles(cycleTable);

%% Plot each gait cycle's time-normalized data, and each gait cycle of one condition plotted on top of each other.
if doPlot
    baseSavePathEMG = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\Filtered_TimeNormalized\EMG';
    plotAllTrials(matchedCycleTable, 'Scaled Time-Normalized EMG', baseSavePathEMG, 'Delsys_Normalized_TimeNormalized');
    baseSavePathXSENS = 'Y:\LabMembers\MTillman\GitRepos\Stroke-R01\Plots\Filtered_TimeNormalized\Joint Angles';
    plotAllTrials(matchedCycleTable, 'Time-Normalized Joint Angles', baseSavePathXSENS, 'XSENS_TimeNormalized');
end

%% Average the data within one SSV/FV & pre/post combination.
avgTableXSENS = avgStructAll(matchedCycleTable, 'XSENS_TimeNormalized', 'XSENS_Averaged', 4);
avgTableDelsys = avgStructAll(matchedCycleTable, 'Delsys_TimeNormalized', 'Delsys_Averaged', 4);
speedPrePostTable = addToTable(speedPrePostTable, avgTableXSENS);
speedPrePostTable = addToTable(speedPrePostTable, avgTableDelsys);

%% SPM Analysis for EMG & XSENS
spmTableXSENS = SPManalysisAll(matchedCycleTable, 'XSENS_TimeNormalized', 'XSENS_SPM', jointsL, jointsR);
spmTableDelsys = SPManalysisAll(matchedCycleTable, 'Delsys_TimeNormalized', 'Delsys_SPM', musclesL, musclesR);
speedPrePostTable = addToTable(speedPrePostTable, spmTableXSENS);
speedPrePostTable = addToTable(speedPrePostTable, spmTableDelsys);

%% Calculate the magnitude and duration of L vs. R differences obtained from SPM in one visit.
magDurTableXSENS = magsDursDiffsLR_All(speedPrePostTable, 'XSENS_SPM', 'XSENS_Averaged', 'XSENS_MagsDiffs');
magDurTableDelsys = magsDursDiffsLR_All(speedPrePostTable, 'Delsys_SPM', 'Delsys_Averaged', 'Delsys_MagsDiffs');
speedPrePostTable = addToTable(speedPrePostTable, magDurTableXSENS);
speedPrePostTable = addToTable(speedPrePostTable, magDurTableDelsys);

%% Calculate root mean squared error (RMSE)
rmseTableXSENS = calculateLRRMSEAll(matchedCycleTable, 'XSENS_TimeNormalized', 'RMSE_JointAngles');
rmseTableDelsys = calculateLRRMSEAll(matchedCycleTable, 'Delsys_Normalized_TimeNormalized', 'RMSE_EMG');
matchedCycleTable = addToTable(matchedCycleTable, rmseTableXSENS);
matchedCycleTable = addToTable(matchedCycleTable, rmseTableDelsys);

%% Calculate cross-correlations
xcorrTableXSENS = calculateLRCrossCorrAll(matchedCycleTable, 'XSENS_TimeNormalized', 'JointAngles_CrossCorr');
xcorrTableDelsys = calculateLRCrossCorrAll(matchedCycleTable, 'Delsys_Normalized_TimeNormalized', 'EMG_CrossCorr');
matchedCycleTable = addToTable(matchedCycleTable, xcorrTableXSENS);
matchedCycleTable = addToTable(matchedCycleTable, xcorrTableDelsys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add this participant's data to the CSV file of all participants' data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Merge the GaitRite and unmatched cycle tables
leadingZeroRegex = '^0+'; % Regex to remove leading zeros
T = copyCategorical(grDistributedTable);
if ismember({'GaitRiteRow'}, T.Properties.VariableNames)
    grDistributedTable.Cycle = grDistributedTable.GaitRiteRow;
    grDistributedTable = removevars(grDistributedTable, 'GaitRiteRow');
end
T = copyCategorical(grDistributedTable);
mergeVarNames = T.Properties.VariableNames(~ismember(T.Properties.VariableNames, {'StartFoot','Cycle'}));

cycleTable.Cycle = categorical(regexprep(string(cycleTable.Cycle),leadingZeroRegex,''));
mergedCycleTable = mergeTables(grDistributedTable, cycleTable, mergeVarNames);

matchedCycleTable.Cycle = categorical(regexprep(string(matchedCycleTable.Cycle),leadingZeroRegex,''));

%% Save the cycle table and the matched cycle table to the all data CSV file
addOneParticipantDataToAllDataCSV(mergedCycleTable, config.PATHS.ALL_DATA_CSV.UNMATCHED);
addOneParticipantDataToAllDataCSV(matchedCycleTable, config.PATHS.ALL_DATA_CSV.MATCHED);
addOneParticipantDataToAllDataCSV(trialTable, config.PATHS.ALL_DATA_CSV.TRIAL);
% addOneParticipantDataToAllDataCSV(visitTable, config.PATHS.ALL_DATA_CSV.VISIT);

%% Save the structs to the participant's save folder.
subjectSavePath = fullfile(subjectSaveFolder, [subject '_' saveFileName]);
if ~isfolder(subjectSaveFolder)
    mkdir(subjectSaveFolder);
end
cycleTable = mergedCycleTable;
save(subjectSavePath, 'visitTable', 'speedPrePostTable', 'trialTable', 'cycleTable', 'matchedCycleTable');
disp(['Saved ' subject ' tables to: ' subjectSavePath]);

