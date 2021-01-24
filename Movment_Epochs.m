%% Movement Epochs
% Examines epochs of movements for rats exploring an open field in between
% stops.

% Stops - epochs of time where the animal moves <3cm/s for at least 2
% seconds


% General maze exploration
save_path = '/Users/lauraberkowitz/Google Drive/Manuscripts/In Progress/TgF344-AD_OF/Data/'; % maxmin.mat, cueCoords.mat

% load data file
load('/Users/lauraberkowitz/Google Drive/Manuscripts/In Progress/TgF344-AD_OF/Data/params.mat'); % maxmin.mat, cueCoords.mat

% Set some measurement parameters
fr = 30; % 30 hz sample rate
diameter = params.dia{1};
binsize=3; %to use for HB detection
upFactor=15; %used in segmentimage - this upscales the heatmap to get a smooth boundary f

% Loop through subjects
for i = 1:length(params.subID)
    
    % Compute measures from back coordinates
    x = params.backCM{i}(:,1);
    y = params.backCM{i}(:,2);
    ts = params.ts{i}(:,1);
    
    % Linear velocity, acceleration, and distance vecotr
    [velocity, acceleration, distance_vector] = linear_motion(x,y,fr,3);
    
    % Get index of when the animal is moving greater than 3cm/second
    move = velocity > 3;
    
    
    %find time when moving greater than or equal to 3cm/s for 1 sec
    runIdx = contiguousframes(move,30);
    [startMotion,endMotion,~] = findgroups(runIdx);
    
    for ii = 1:length(startMotion)
        motion{ii} = [x(startMotion(ii):endMotion(ii)),y(startMotion(ii):endMotion(ii))];
        timeMotion{ii} = size(motion{ii},1)/fr;
        tsMotion{ii} = ts(startMotion(ii):endMotion(ii),1);
    end
    
    runs = motion;   clear motion
    time_moving = timeMotion; clear timeMotion
    ts_motion  = tsMotion; clear tsMotion
    num_runs  = size(runs ,2);
        
    %Compute circuity between segments
    for iii = 1:size(runs,2)
        epoch_circuity(iii) = circuity(runs{iii}(:,1),runs{iii}(:,2));
    end
    
    %
    mean_circ = nanmedian(epoch_circuity(1,:)); %compute average circuity
    num_runs = size(runs,2);
    
%     stop_measures = OF.stops(x,y,ts,velocity,fr,2*fr);

    % Subject, group, and condition identifiers
    subID{i,1} = params.subID{i};
    temp_group = contains(subID{i,1},'Tg');
    if temp_group == 0
        group{i,1} = 'wt';
    else
        group{i,1} = 'tg';
    end
    day{i,1} = extractAfter(params.subID{i},'_');
    
end

% Save to google drive
vars = {'subID','group','day','path_length',...
    'median_velocity',...
    'median_abs_acceleration',...
    'median_abs_angular_velocity',...
    'search_area'};

df = table(subID,group,day,path_length,median_velocity,median_acceleration,...
    median_angular_velocity,search_area,'VariableNames',vars);


% Writetable
writetable(df,[save_path,'general_locomotion.csv'])

