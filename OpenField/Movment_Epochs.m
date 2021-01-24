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
    cue_cords = params.cueCM{i};  
    
    % Linear velocity, acceleration, and distance vecotr
    [velocity, acceleration, distance_vector] = linear_motion(x,y,fr,3);
    
    % Get index of when the animal is moving greater than 3cm/second
    move = velocity > 3;
    
    
    %find time when moving greater than or equal to 3cm/s for 1 sec
    runIdx = contiguousframes(move,30);
    [startMotion,endMotion,~] = findgroups(runIdx);
    
    for ii = 1:length(startMotion)
        runs{ii} = [x(startMotion(ii):endMotion(ii)),y(startMotion(ii):endMotion(ii))];
        time_moving{ii} = size(runs{ii},1)/fr;
        ts_motion{ii} = ts(startMotion(ii):endMotion(ii),1);
    end
    
        
    %Compute circuity and cue proximity, between run segments
    for iii = 1:size(runs,2)
        run_circuity(iii) = circuity(runs{iii}(:,1),runs{iii}(:,2));
        run_cue_proximity(iii)  = proximity_to_cue(runs{iii}(:,1),runs{iii}(:,2),cue_cords);
        run_velocity(iii) = nanmedian(velocity(ismember(ts,ts_motion{iii})));
    end
    
    % 
    circ = nanmedian(run_circuity); %compute median circuity
    num_runs = size(runs,2);
    proximity_cue = nanmedian(run_cue_proximity); 
    
    
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


%% Local functions

% Movement proximity to cue

function move_cue_proximity = proximity_to_cue(x,y,cue_cords)
cue_x = cue_cords(:,1);
cue_y = cue_cords(:,2);

% convex hull to close shape of cue
k = convhull(cue_x,cue_y);
cue_boundary = [cue_x(k) cue_y(k)];

% distance from cue boundary
for r = 1:length(x)
    distances = sqrt(sum(bsxfun(@minus, cue_boundary, [x(r),y(r)]).^2,2));
    move_cue_distance(r) = unique(distances(distances==min(distances))); %Find minimum distance from cue boundary to hb center
end

% find median distance 
move_cue_proximity = median(move_cue_distance);

end

