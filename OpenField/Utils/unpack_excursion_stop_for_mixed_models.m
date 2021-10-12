%% Unpack params for to long format for linear mixed effects modelling 

% General maze exploration
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/'; % maxmin.mat, cueCoords.mat
% load data file
load('/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/params.mat'); % maxmin.mat, cueCoords.mat

%%
day = []; group = []; subID = []; runID = []; run_circ = []; run_duration = []; run_length = []; cue_bearing = [];
run_peak_vel = []; run_cue_proximity = []; stop_time = []; stop_cue_proximity = [];
for i = 1:length(params.subID)
    
    % repmat day, group, subID for length of run segments 
    sub_group = contains(params.subID{i,1},'Tg');
    if sub_group == 0
        temp_group = 'wt';
    else
        temp_group = 'tg';
    end
    to_rep = size(params.runs{i},2);
    day = [day; repmat(extractAfter(params.subID{i},'_'),to_rep,1)];
    group = [group; repmat(temp_group,to_rep,1)];
    subID = [subID; repmat(extractBefore(params.subID{i},'_'),to_rep,1)];
    runID = [runID; [1:to_rep]'];
    
    % concatenate data 
    run_circ = [run_circ; params.run_circ{i}'];
    run_duration = [run_duration; params.run_duration{i}];
    run_length = [run_length; params.run_length{i}'];
    cue_bearing = [cue_bearing; params.cue_bearing{i}'];
    run_peak_vel = [run_peak_vel; params.peak_vel{i}];
    run_cue_proximity = [run_cue_proximity; params.cue_proximity{i}'];
    
end

%%

% Save to google drive
vars = {'subID','group','day','run_id','circuity',...
    'duration',...
    'length',...
    'cue_bearing',...
    'run_peak_vel',...
    'run_cue_proximity',...
    };

df = table(subID,group,day,runID,run_circ,run_duration,run_length,cue_bearing,run_peak_vel,run_cue_proximity,'VariableNames',vars);


% Writetable
writetable(df,[save_path,'movement_metrics_mixed_models.csv'])

%% Stop analysis 
day = []; group = []; subID = []; stopID = [];
stop_time = []; stop_cue_proximity = []; inter_stop_interval =[]; stop_start_time = [];
for i = 1:length(params.subID)
    
    % repmat day, group, subID for length of run segments 
    sub_group = contains(params.subID{i,1},'Tg');
    if sub_group == 0
        temp_group = 'wt';
    else
        temp_group = 'tg';
    end
    to_rep = size(params.timeStopped{i},1);
    day = [day; repmat(extractAfter(params.subID{i},'_'),to_rep,1)];
    group = [group; repmat(temp_group,to_rep,1)];
    subID = [subID; repmat(extractBefore(params.subID{i},'_'),to_rep,1)];
    stopID = [stopID; [1:to_rep]'];
    
    % concatenate data 
    stop_time = [stop_time; params.timeStopped{i}];
    stop_cue_proximity = [stop_cue_proximity; params.stop_cue_proximity{i}'];
    inter_stop_interval = [inter_stop_interval; params.inter_stop_interval{i}];
    stop_start_time = [stop_start_time; params.stop_start{i}'];

end

vars = {'subID','group','day','stop_id','stop_time','stop_cue_proximity','inter_stop_interval','stop_start_time'};

df = table(subID,group,day,stopID,stop_time,stop_cue_proximity,inter_stop_interval,stop_start_time,'VariableNames',vars);

% Writetable
writetable(df,[save_path,'stop_metrics_mixed_models.csv'])

