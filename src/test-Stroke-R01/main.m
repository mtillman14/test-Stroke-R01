paths; % Get the paths
config = jsondecode(fileread(configPath));

runConfig = toml.map_to_struct(toml.read(subjectsToRunPath));
allSubjects = runConfig.subjects.run;

%% Iterate over each subject
doPlot = false;
for subNum = 1:length(allSubjects)
    subject = allSubjects{subNum};    
    subjectSavePath = fullfile(config.PATHS.ROOT_SAVE, subject, [subject '_' config.PATHS.SAVE_FILE_NAME]);
    disp(['Now running subject (' num2str(subNum) '/' num2str(length(allSubjects)) '): ' subject]);
    mainOneSubject; % Run the main pipeline.
end

%% Plot each subject
allSubjectsPlot = runConfig.subjects.plot;
for subNum = 1:length(allSubjectsPlot)
    subject = allSubjectsPlot{subNum};
    loadPath = fullfile(config.PATHS.ROOT_SAVE, subject, [subject '_Overground_EMG_Kinematics.mat']);
    load(loadPath, 'matchedCycleTable');
    % Plot each gait cycle's filtered data, time normalized (for EMG, scaled to max EMG) and each gait cycle of one condition plotted on top of each other.
    baseSavePath = fullfile(config.PATHS.PLOTS.ROOT, config.PATHS.PLOTS.FILTERED_TIME_NORMALIZED);
    baseSavePathEMG = fullfile(baseSavePath, 'EMG');
    baseSavePathXSENS = fullfile(baseSavePath, 'Joint Angles');
    % plotAllTrials(matchedCycleTable, 'Time-Normalized Non-Normalized EMG', baseSavePathEMG, 'Delsys_TimeNormalized'); 
    plotAllTrials(matchedCycleTable, 'Time-Normalized Scaled EMG', baseSavePathEMG, 'Delsys_Normalized_TimeNormalized'); 
    % plotAllTrials(matchedCycleTable, 'Time-Normalized Joint Angles', baseSavePathXSENS, 'XSENS_TimeNormalized');
end

%% Load the cycleTable and matchedCycleTable from all subjects
paths;
config = jsondecode(fileread(configPath));
categoricalCols = {'Subject','Intervention','PrePost','Speed','Trial','Cycle','StartFoot'};
cycleTableAll = readtable(config.PATHS.ALL_DATA_CSV.UNMATCHED);
matchedCycleTableAll = readtable(config.PATHS.ALL_DATA_CSV.MATCHED);
for i = 1:length(categoricalCols)
    cycleTableAll.(categoricalCols{i}) = categorical(cycleTableAll.(categoricalCols{i}));
    matchedCycleTableAll.(categoricalCols{i}) = categorical(matchedCycleTableAll.(categoricalCols{i}));
end

%% Load the trialTable from all subjects
trialTableCategoricalCols = {'Subject','Intervention','PrePost','Speed','Trial'};
trialTableAll = readtable(config.PATHS.ALL_DATA_CSV.TRIAL);
for i = 1:length(trialTableCategoricalCols)
    trialTableAll.(trialTableCategoricalCols{i}) = categorical(trialTableAll.(trialTableCategoricalCols{i}));
end

%% Replace 'StartFoot' with 'Side'
categoricalCols = {'Subject','Intervention','PrePost','Speed','Trial','Cycle','Side'};
if ismember('StartFoot', cycleTableAll.Properties.VariableNames)
    cycleTableAll.Side = cycleTableAll.StartFoot;
    cycleTableAll = removevars(cycleTableAll, 'StartFoot');
end
if ismember('StartFoot', matchedCycleTableAll.Properties.VariableNames)
    matchedCycleTableAll.Side = matchedCycleTableAll.StartFoot;
    matchedCycleTableAll = removevars(matchedCycleTableAll, 'StartFoot');
end
cycleTableAll = movevars(cycleTableAll,'Side','After','Cycle');
matchedCycleTableAll = movevars(matchedCycleTableAll,'Side','After','Cycle');

%% Calculate symmetries
formulaNum = 6; % computing Sym using the original equation * 100
[colNamesL, colNamesR] = getLRColNames(cycleTableAll);
% Cycle table
cycleTableContraRemoved_NoGR = removeContralateralSideColumns(cycleTableAll, colNamesL, colNamesR);
grVars = cycleTableAll.Properties.VariableNames(contains(cycleTableAll.Properties.VariableNames,'_GR'));
grTable = removevars(cycleTableAll, ~ismember(cycleTableAll.Properties.VariableNames, [grVars, categoricalCols]));
cycleTableContraRemoved = addToTable(cycleTableContraRemoved_NoGR, grTable);
scalarColumnNames = getScalarColumnNames(cycleTableContraRemoved);
allColumnNames = cycleTableContraRemoved.Properties.VariableNames;
nonscalarColumnNames = allColumnNames(~ismember(allColumnNames, [scalarColumnNames; categoricalCols']));
cycleTableContraRemovedScalarColumns = removevars(cycleTableContraRemoved, nonscalarColumnNames);
% Compute the symmetry values
nonSubsetCatVars = {'Cycle','Side'};
lrSidesCycleSymTable = calculateSymmetryAll(cycleTableContraRemovedScalarColumns, '_Sym', formulaNum, nonSubsetCatVars);
categoricalColsTrial = {'Subject','Intervention','PrePost','Speed','Trial'};
trialTableAllSym = trialTableAll;
cycleTableAllSym = cycleTableContraRemovedScalarColumns;
matchedCycleTableAllSym = addToTable(matchedCycleTableAll, lrSidesCycleSymTable);

%% Adjust intervention name to mapped names
mapped_interventions = config.MAPPED_INTERVENTION_FIELDS;
interventions = config.INTERVENTION_FOLDERS;
intervention_map = containers.Map(interventions, mapped_interventions);
for i = 1:height(trialTableAllSym)
    trialTableAllSym.Intervention(i) = intervention_map(string(trialTableAllSym.Intervention(i)));
end
for i = 1:height(cycleTableAllSym)
    cycleTableAllSym.Intervention(i) = intervention_map(string(cycleTableAllSym.Intervention(i)));
end
for i = 1:height(matchedCycleTableAllSym)
    matchedCycleTableAllSym.Intervention(i) = intervention_map(string(matchedCycleTableAllSym.Intervention(i)));
end

%% Add the StimNoStim, Intensity, and Frequency columns
interventionColumnName = 'Intervention';
trialTableAllAddedCols = addStimNoStim_Intensity_FrequencyCols(trialTableAllSym, interventionColumnName);
cycleTableAllAddedCols = addStimNoStim_Intensity_FrequencyCols(cycleTableAllSym, interventionColumnName);
matchedCycleTableAllAddedCols = addStimNoStim_Intensity_FrequencyCols(matchedCycleTableAllSym, interventionColumnName);

%% Add session number
addpath('Y:\LabMembers\MTillman\GitRepos\Stroke-R01\src\MEPs\MEPs Processing AIM 1');
tepsLogPath = 'Y:\Spinal Stim_Stroke R01\AIM 1\Subject Data\TEPs_log.xlsx';
tepsLog = readExcelFileOneSheet(tepsLogPath, 'Subject','Sheet1');
allColNames = tepsLog.Properties.VariableNames;
colNames = {'Subject', 'SessionOrder', 'SessionCode'};
colNamesIdx = ismember(allColNames, colNames);
reducedTEPsLog = unique(tepsLog(:, colNamesIdx), 'rows');
for i = 1:height(reducedTEPsLog)
    reducedTEPsLog.Subject{i} = ['SS' reducedTEPsLog.Subject{i}];
end
subjectIdx = ismember(string(reducedTEPsLog.Subject), string(trialTableAll.Subject));
reducedTEPsLog(~subjectIdx,:) = [];
% Map the intervention names
mappedInterventions = containers.Map(config.INTERVENTION_FOLDERS, config.MAPPED_INTERVENTION_FIELDS);
reducedTEPsLog.SessionCode = cellfun(@(x) mappedInterventions(x), reducedTEPsLog.SessionCode, 'UniformOutput', false);
sessionOrderColName = 'SessionOrder';
sessionCodeColName = 'SessionCode';
interventionColName = 'Intervention';
trialTableAllSessionNum = addSessionOrder(trialTableAllAddedCols, reducedTEPsLog, sessionOrderColName, sessionCodeColName, interventionColName, interventionColName);
cycleTableAllSessionNum = addSessionOrder(cycleTableAllAddedCols, reducedTEPsLog, sessionOrderColName, sessionCodeColName, interventionColName, interventionColName);
matchedCycleTableAllSessionNum = addSessionOrder(matchedCycleTableAllAddedCols, reducedTEPsLog, sessionOrderColName, sessionCodeColName, interventionColName, interventionColName);

%% Adjust the L & R sides to "U" and "A" for unaffected and affected sides
tepsLogPath = 'Y:\Spinal Stim_Stroke R01\AIM 1\Subject Data\TEPs_log.xlsx';
tepsLog = readExcelFileOneSheet(tepsLogPath, 'Subject','Sheet1');
colNames = {'Subject','PareticSide'};
inputTableSideCol = 'Side';
tepsLogSideCol = 'PareticSide';
allColNames = tepsLog.Properties.VariableNames;
colNamesIdx = ismember(allColNames, colNames);
reducedTEPsLog = unique(tepsLog(:, colNamesIdx), 'rows');
for i = 1:height(reducedTEPsLog)
    reducedTEPsLog.Subject{i} = ['SS' reducedTEPsLog.Subject{i}];
end
subjectIdx = ismember(string(reducedTEPsLog.Subject), string(trialTableAll.Subject));
reducedTEPsLog(~subjectIdx,:) = [];
trialTableAllUA = trialTableAllSessionNum;
cycleTableAllUA = convertLeftRightSideToAffectedUnaffected(cycleTableAllSessionNum, reducedTEPsLog, inputTableSideCol, tepsLogSideCol);
matchedCycleTableAllUA = convertLeftRightSideToAffectedUnaffected(matchedCycleTableAllSessionNum, reducedTEPsLog, inputTableSideCol, tepsLogSideCol);

%% Save the unaffected and affected side tables
tablesPathPrefixMergedUA = 'Y:\LabMembers\MTillman\SavedOutcomes\StrokeSpinalStim\Overground_EMG_Kinematics\MergedTablesAffectedUnaffected_0_1_reproduced';
writetable(trialTableAllUA, fullfile(tablesPathPrefixMergedUA, 'trialTableAll.csv'));
writetable(matchedCycleTableAllUA, fullfile(tablesPathPrefixMergedUA, 'matchedCycles.csv'));
writetable(cycleTableAllUA, fullfile(tablesPathPrefixMergedUA, 'unmatchedCycles.csv'));

%% Calculate pre to post change
categoricalColsNoPrePost = {'Subject','Intervention','Speed','Trial','Cycle','Side'};
% Percent difference
formulaNum = 2;
prePostCycleChangeTablePercDiff = calculatePrePostChange(cycleTableAllUA, formulaNum);
prePostChangeMatchedCycleTablePercDiff = calculatePrePostChange(matchedCycleTableAllUA, formulaNum);
% Difference
formulaNum = 1;
prePostCycleChangeTableDiff = calculatePrePostChange(cycleTableAllUA, formulaNum);
prePostChangeMatchedCycleTableDiff = calculatePrePostChange(matchedCycleTableAllUA, formulaNum);
% Combine the two tables 
prePostCycleChangeTable = join(prePostCycleChangeTableDiff, prePostCycleChangeTablePercDiff, 'Keys', categoricalColsNoPrePost);
prePostChangeMatchedCycleTable = join(prePostChangeMatchedCycleTableDiff, prePostChangeMatchedCycleTablePercDiff, 'Keys', categoricalColsNoPrePost);

%% Add the 10MWT data to each table
trialTableAllSessionNum10MWT = join10MWTSpeedToCycleLevelTable(tepsLogPath, fullfile(tablesPathPrefixMergedUA, 'trialTableAll.csv'), configPath);
mergedMatchedCycleTableUA10MWT = join10MWTSpeedToCycleLevelTable(tepsLogPath, fullfile(tablesPathPrefixMergedUA, 'matchedCycles.csv'), configPath);
mergedUnmatchedCycleTableUA10MWT = join10MWTSpeedToCycleLevelTable(tepsLogPath, fullfile(tablesPathPrefixMergedUA, 'unmatchedCycles.csv'), configPath);

%% Save the 10MWT tables
tablesPathPrefixMergedUA10MWT = 'Y:\LabMembers\MTillman\SavedOutcomes\StrokeSpinalStim\Overground_EMG_Kinematics\MergedTablesAffectedUnaffected10MWT';
writetable(trialTableAllSessionNum10MWT, fullfile(tablesPathPrefixMergedUA10MWT, 'trialTableAll.csv'));
writetable(mergedMatchedCycleTableUA10MWT, fullfile(tablesPathPrefixMergedUA10MWT, 'matchedCycles.csv'));
writetable(mergedUnmatchedCycleTableUA10MWT, fullfile(tablesPathPrefixMergedUA10MWT, 'unmatchedCycles.csv'));

%% Widen the unaffected and affected side matchedCycle tables
inputTableSideCol = 'Side';
factorColNames = {'Subject','Intervention','Speed','Trial', 'PrePost'};
preStruct.PrePost = 'PRE';
postStruct.PrePost = 'POST';
mergedMatchedCycleTablePath = 'Y:\LabMembers\MTillman\SavedOutcomes\StrokeSpinalStim\Overground_EMG_Kinematics\MergedTablesAffectedUnaffected10MWT\matchedCycles.csv';
mergedMatchedCycleTable = readtable(mergedMatchedCycleTablePath);
mergedMatchedCycleTableUAWidePreMean = widenTableBySides(mergedMatchedCycleTable, inputTableSideCol, factorColNames, preStruct, 'mean');
mergedMatchedCycleTableUAWidePreMedian = widenTableBySides(mergedMatchedCycleTable, inputTableSideCol, factorColNames, preStruct, 'median');
mergedMatchedCycleTableUAWidePostMean = widenTableBySides(mergedMatchedCycleTable, inputTableSideCol, factorColNames, postStruct, 'mean');
mergedMatchedCycleTableUAWidePostMedian = widenTableBySides(mergedMatchedCycleTable, inputTableSideCol, factorColNames, postStruct, 'median');

%% Save the widened matchedCycle tables
mergedMatchedCycleTableWidePathPrefix = 'Y:\LabMembers\MTillman\SavedOutcomes\StrokeSpinalStim\Overground_EMG_Kinematics\MergedTablesAffectedUnaffectedWide';
writetable(mergedMatchedCycleTableUAWidePreMean, fullfile(mergedMatchedCycleTableWidePathPrefix, 'matchedCycles_pre_mean.csv'));
writetable(mergedMatchedCycleTableUAWidePreMedian, fullfile(mergedMatchedCycleTableWidePathPrefix, 'matchedCycles_pre_median.csv'));
writetable(mergedMatchedCycleTableUAWidePostMean, fullfile(mergedMatchedCycleTableWidePathPrefix, 'matchedCycles_post_mean.csv'));
writetable(mergedMatchedCycleTableUAWidePostMedian, fullfile(mergedMatchedCycleTableWidePathPrefix, 'matchedCycles_post_median.csv'));

%% Widen the unaffected and affected side unmatchedCycle tables
mergedUnmatchedCycleTablePath = 'Y:\LabMembers\MTillman\SavedOutcomes\StrokeSpinalStim\Overground_EMG_Kinematics\MergedTablesAffectedUnaffected10MWT\unmatchedCycles.csv';
mergedUnmatchedCycleTable = readtable(mergedUnmatchedCycleTablePath);
mergedUnmatchedCycleTableUAWidePreMean = widenTableBySides(mergedUnmatchedCycleTable, inputTableSideCol, factorColNames, preStruct, 'mean');
mergedUnmatchedCycleTableUAWidePreMedian = widenTableBySides(mergedUnmatchedCycleTable, inputTableSideCol, factorColNames, preStruct, 'median');
mergedUnmatchedCycleTableUAWidePostMean = widenTableBySides(mergedUnmatchedCycleTable, inputTableSideCol, factorColNames, postStruct, 'mean');
mergedUnmatchedCycleTableUAWidePostMedian = widenTableBySides(mergedUnmatchedCycleTable, inputTableSideCol, factorColNames, postStruct, 'median');

%% Save the widened unmatchedCycle tables
mergedMatchedCycleTableWidePathPrefix = 'Y:\LabMembers\MTillman\SavedOutcomes\StrokeSpinalStim\Overground_EMG_Kinematics\MergedTablesAffectedUnaffectedWide';
writetable(mergedUnmatchedCycleTableUAWidePreMean, fullfile(mergedMatchedCycleTableWidePathPrefix, 'unmatchedCycles_pre_mean.csv'));
writetable(mergedUnmatchedCycleTableUAWidePreMedian, fullfile(mergedMatchedCycleTableWidePathPrefix, 'unmatchedCycles_pre_median.csv'));
writetable(mergedUnmatchedCycleTableUAWidePostMean, fullfile(mergedMatchedCycleTableWidePathPrefix, 'unmatchedCycles_post_mean.csv'));
writetable(mergedUnmatchedCycleTableUAWidePostMedian, fullfile(mergedMatchedCycleTableWidePathPrefix, 'unmatchedCycles_post_median.csv'));