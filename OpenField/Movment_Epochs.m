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
binsize = 3; %to use for HB detection
upFactor = 15; %used in segmentimage - this upscales the heatmap to get a smooth boundary f

% Save variables 
params.runs{1} = [];
params.run_circ{1} = [];
% Loop through subjects
for i = 1:length(params.subID)
    
    % Compute measures from back coordinates
    x = params.backCM{i}(:,1);
    y = params.backCM{i}(:,2);
    ts = params.ts{i}(:,1);
    cue_cords = params.cueCM{i};
    
    % Linear velocity, acceleration, and distance vecotr
    [velocity, acceleration, distance_vector] = linear_motion(x,y,fr,0.25);
    vel_ts = 0+(1/fr)/2:1/fr:max(ts)-(1/fr)/2;
    
    idx = isnan(velocity);
    velocity(idx) = [];
    vel_ts(idx) = [];
    
    idx = isnan(x) | isnan(y);
    x(idx) = [];
    y(idx) = [];
    ts(idx) = [];
%     velocity((isnan(velocity))) = [];
%     vel_ts((isnan(velocity))) = [];
    %%
    % Find events when rat is running
    [events,duration,peak_vel] = detect_velocity_epochs(velocity,vel_ts,fr);
    %%
    % Grab xy for each run
    for ii = 1:length(events)
        idx = ts >= events(ii,1) & ts <= events(ii,2);
        run{ii} = [x(idx),y(idx)];
    end
    
    
    %Compute circuity and cue proximity, between run segments
    for iii = 1:size(run,2)
        run_circuity(iii) = circuity(run{iii}(:,1),run{iii}(:,2));
        run_length(iii) = nansum(sqrt((diff(run{iii}(:,1))).^2 + (diff(run{iii}(:,2))).^2));

        if isnan(cue_cords)
            run_cue_proximity(iii)  = proximity_to_cue(run{iii}(:,1),run{iii}(:,2),params.cueCM{i - 1});
        else
            run_cue_proximity(iii)  = proximity_to_cue(run{iii}(:,1),run{iii}(:,2),cue_cords);
            
        end
    end
    
    % save run coords, circuity, and run length back to params for easy
    % plotting. 
    params.runs{i} = run;
    params.run_circ{i} = run_circuity;
    params.run_length{i} = run_length;
    clear run
    
    % Calulate some gross run metrics
    circ(i) = nanmedian(run_circuity); %compute median circuity
    num_runs(i)  = size(duration,1);
    proximity_cue(i)  = nanmedian(run_cue_proximity);
    median_peak_vel(i)  = nanmedian(peak_vel);
    median_duration(i)  = nanmedian(duration);
    median_inter_run_interval(i)  = median(diff(events(:,3)));
    median_run_length(i) = median(run_length);
    
    
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
vars = {'subID','group','day','circuity',...
    'num_runs',...
    'run_length',...
    'proximity_cue',...
    'median_peak_vel',...
    'median_duration',...
    'median_inter_run_interval',...
    };

df = table(subID,group,day,circ',num_runs',median_run_length',proximity_cue',...
    median_peak_vel',median_duration',median_inter_run_interval','VariableNames',vars);


% Writetable
writetable(df,[save_path,'movement_epochs.csv'])


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

function [events,duration,peak_vel] = detect_velocity_epochs(velocity,vel_ts,fr)
% RH
velocity_thres_sd = 0; % standard deviations
highThresholdFactor = 2;
min_inter_event_interval = 1; % seconds
min_event_duration = 1; % seconds

velocity = zscore(velocity);

thresholded = velocity > velocity_thres_sd;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);

% Exclude last epoch if it is incomplete
if length(stop) == length(start)-1
    start = start(1:end-1);
end
% Exclude first epoch if it is incomplete
if length(stop)-1 == length(start)
    stop = stop(2:end);
end
% Correct special case when both first and last epoch are incomplete
if start(1) > stop(1)
    stop(1) = [];
    start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass)
    disp('Detection by thresholding failed');
    return
else
    disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

% Merge event if inter-event period is too short
min_inter_event_samples = min_inter_event_interval*fr;
secondPass = [];
event = firstPass(1,:);
for i = 2:size(firstPass,1)
    if firstPass(i,1) - event(2) < min_inter_event_samples
        % Merge
        event = [event(1) firstPass(i,2)];
    else
        secondPass = [secondPass ; event];
        event = firstPass(i,:);
    end
end
secondPass = [secondPass ; event];
if isempty(secondPass)
    disp('event merge failed');
    return
else
    disp(['After event merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard event with a peak power < highThresholdFactor
thirdPass = [];
peak_vel = [];
for i = 1:size(secondPass,1)
    [maxValue,maxIndex] = max(velocity([secondPass(i,1):secondPass(i,2)]));
    if maxValue > highThresholdFactor
        thirdPass = [thirdPass ; secondPass(i,:)];
        peak_vel = [peak_vel ; maxValue];
    end
end
if isempty(thirdPass)
    disp('Peak thresholding failed.');
    return
else
    disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end

% Detect peak position for each event
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1)
    [minValue,minIndex] = max(velocity(thirdPass(i,1):thirdPass(i,2)));
    peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end

events = [vel_ts(thirdPass(:,1))' ,vel_ts(thirdPass(:,2))',vel_ts(peakPosition)'];
duration = events(:,2)-events(:,1);

% Discard events that are too short
events(duration < min_event_duration,:) = NaN;
events = events((all((~isnan(events)),2)),:);
duration(duration < min_event_duration,:) = NaN;
duration = duration((all((~isnan(duration)),2)),:);

disp(['After duration test: ' num2str(size(events,1)) ' events.']);
%%
%     figure;
%     plot(vel_ts,velocity)
%     hold on
%     for ep = 1:length(events)
%         plot([events(ep,1),events(ep,1)],ylim,'g','LineWidth',3)
%         plot([events(ep,2),events(ep,2)],ylim,'r','LineWidth',3)
%     end
%
%     figure
%     histogram(diff(events(:,3)),200)

end

function length = path_length(x,y,velocity)
% Computes total path length during running
% inputs:
%   x: position vector of x-coordinates of length n
%   y: position vecotr of y-coordinates of length n
%   velocity: vector of instantaneous velocity (length of n - 1)
% output:
%   length: total path length in cm for active motion


%distance formula
distance_vector = sqrt((diff(x)).^2 + (diff(y)).^2);

% Summary Path Measures
length = sum(distance_vector(velocity >= 3,1));%Total Path length for points greater than 3cm/s

end

