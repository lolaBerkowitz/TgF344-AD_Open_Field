% small_openfield_analysis examines locomotion and home base behaviors of
% TgF344-AD and F344 controls in small open field (76.5cm) in 10min
% sessions across days (days 1,2,3) and time-points (4 months, 7 months, 10
% months). Positional data was tracked using anymaze and xy coordinates
% were saved as .csv files using the Process_SLK.ipynb

% Save path (where final csv is saved)
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/';

% video acquisition parameters
fr = 9; % frame rate in Hz
diameter = 76.5; % diameter of circular maze
binsize=3; %to use for HB detection
upFactor=15; %used in segmentimage - this upscales the heatmap to get a smooth boundary f

% initalize data table
df = table;

% Import CSV files
files = dir("/Users/lauraberkowitz/Desktop/coords_csv/**/*.csv");

% Pull identifiers from file name
for i = 1:length(files)
    temp = strsplit(files(i).name,' ');
    df.file_loc{i} = [files(i).folder,filesep,files(i).name];
    df.genotype{i} = temp{1};
    df.rat(i) = str2double(temp(2));
    df.day(i) = str2double(temp(4));
    df.time_point{i} = sscanf(temp{end},'%d');
end

clear temp files
%% Preprocess paths (estimate max/min critical for transform into cm)
% add estimated maxmin to df
df = estimate_maxmin(df);

% transform into cm and smooths paths
df = transform_all_coords(df,diameter);


%% Locomotion (path length, median running speed, occupancy, thigmotaxis)
df = locomotion_metrics(df,fr,diameter,binsize);

%% Stops (stop frequency, mean/median stop duration)
df = stop_metrics(df,fr);

%% Homebase (number of home bases)
df = home_base_analysis(df,fr,binsize,diameter,15);

% plot distributions for duration by area
% fig = figure; 
% fig.Color = [1 1 1];
% subplot(1,2,1)
% histogram([df.hb_duration_all{:}],100); hold on;
% xline(50,'r','LineWidth',2)
% xlabel('Home base duration (sec)')
% ylabel('probability')
% 
% subplot(1,2,2)
% scatter([df.hb_area{:}],[df.hb_duration_all{:}]); hold on;
% xlabel('Home base area (sq cm)')
% ylabel('Home base duration (sec)')
% 
% for hb = 1:length(df.hb_duration_all)
%     
%     [temp_duration,idx] = max(df.hb_duration_all{hb});
%     temp_area = df.hb_area{hb}(idx);
%     
%     scatter(temp_area,temp_duration,'filled','r')
%    
% end
% legend({'home base','primary home base'})
% 
% export_fig('/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/figs/small_hb_testing.png','-m4') 

%% Home base distance by timepoint 
load('df_small_OF.mat');
primary_hb_dist_mat = hb_dist_by_day(df);

vars = {'rat','time_point','group','mean_hb_distance_over_days'};

for i = 1:length(primary_hb_dist_mat)
    if primary_hb_dist_mat(i,1) < 200 
        group{i,1} = 'tg';
    else
        group{i,1} = 'wt';
    end
end

primary_hb_dist_table = table(primary_hb_dist_mat(:,1),primary_hb_dist_mat(:,2),...
    group,primary_hb_dist_mat(:,3),'VariableNames',vars);

writetable(primary_hb_dist_table,[save_path,filesep,'primary_hb_dist_table.csv'])


%% Save outcome measures to csv for analysis in R
unpack_df(df,save_path)

%% save cm paths
for i = 1:length(df.file_loc)
    temp_table = table;
    temp_table.ts = df.ts{i};
    temp_table.x = df.x_cm{i};
    temp_table.y = df.y_cm{i};
    temp_vel = df.velocity{i};
    temp_vel(temp_vel < 0) = 0;
    temp_table.velocity = temp_vel;
    temp_id = [df.genotype{i},'_',num2str(df.rat(i)),'_',num2str(df.time_point{i}),'_',num2str(df.day(i))];
    writetable(temp_table,['/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/small_of_paths/',temp_id,'.csv'])
end

%% local functions

function unpack_df(df,save_path)
%% Stop analysis 
day = []; group = []; subID = []; time_point = [];

% locomotion measures
gen_loco_thigmotaxis = []; gen_loco_velocity = []; gen_loco_path_length = []; gen_loco_search_area = [];
% stop measures
stops_num_stops = []; stops_median_duration = []; stops_median_inter_stop_interval = []; 
% home base
hb_duration = []; hb_stop_dist = []; hb_close_stop = []; hb_time2HB_keep = []; hb_HB_count_keep = []; hb_dist_to_wall = []; 

% homebase measures
for i = 1:length(df.file_loc)
    
    % repmat day, group, subID for length of run segments 
    sub_group = contains(df.genotype{i,1},'Tg');
    if sub_group == 0
        temp_group = 'wt';
    else
        temp_group = 'tg';
    end
    
    % factor identifiers 
    day = [day; df.day(i)];
    time_point = [time_point; df.time_point(i)];
    group = [group; temp_group];
    subID = [subID; df.rat(i)];
    
    % concatenate data 
    gen_loco_thigmotaxis = [gen_loco_thigmotaxis; df.gen_loco_thigmotaxis(i)];
    gen_loco_velocity = [gen_loco_velocity; df.gen_loco_velocity(i)];
    gen_loco_path_length = [gen_loco_path_length; df.gen_loco_path_length(i)];
    gen_loco_search_area = [gen_loco_search_area; df.gen_loco_search_area(i)];

    stops_num_stops = [stops_num_stops; df.stops_num_stops(i)];
    stops_median_duration = [stops_median_duration; df.stops_median_duration(i)];
    stops_median_inter_stop_interval = [stops_median_inter_stop_interval; df.stops_median_inter_stop_interval(i)];

    hb_duration = [hb_duration; df.hb_duration(i)];
    hb_stop_dist = [hb_stop_dist; df.hb_stop_dist(i)];
    hb_close_stop = [hb_close_stop; df.hb_close_stop(i)];
    hb_time2HB_keep = [hb_time2HB_keep; df.hb_time2HB_keep(i)];
    hb_HB_count_keep = [hb_HB_count_keep; df.hb_HB_count_keep(i)];
    hb_dist_to_wall = [hb_dist_to_wall; df.hb_min_dist_from_wall(i)];


end

vars = {'subID','group','day','time_point',...
    'thigmotaxis','median_velocity','path_length','search_area'...
    'n_stops','median_stop_duration','median_inter_stop_interval'...
    'time_in_homebase','hb_stop_distance','hb_close_stops','time2hb','n_home_bases','hb_dist_to_wall'};

save_df = table(subID,group,day,time_point,...
    gen_loco_thigmotaxis,gen_loco_velocity,gen_loco_path_length,gen_loco_search_area,...
    stops_num_stops,stops_median_duration,stops_median_inter_stop_interval,...
    hb_duration,hb_stop_dist,hb_close_stop,hb_time2HB_keep,hb_HB_count_keep,hb_dist_to_wall,'VariableNames',vars);

% Writetable
writetable(save_df,[save_path,filesep,'small_open_field_measures.csv'])
disp(['data compiled to csv and saved at :',save_path])
end

% quick and dirty way to get max min from paths
function [x_max,x_min,y_max,y_min] = path_maxmin(x,y)
x_max = max(x);
x_min = min(x);
y_max = max(y);
y_min = min(y);
end

function df = estimate_maxmin(df)
% estimates the boundary of the apparatus by calculating the max/min of
% coordinates for all paths within a time point. depends on at least one
% animal exploring the edges of an apparatus & the apparatus not moving
% between trials.
df.raw_x{1} = [];
df.raw_y{1} = [];
%loop through df files and pull max/min
for i = 1:length(df.file_loc)
    temp = readtable(df.file_loc{i});
    x = temp.x;
    y = temp.y;
    ts = temp.ts;
    % save coordinates back to data table for ease
    df.raw_x{i} = x;
    df.raw_y{i} = y;
    df.ts{i} = ts;
    [x_max(i,1),x_min(i,1),y_max(i,1),y_min(i,1)] = path_maxmin(x,y);
    
    
end

% for given time point, find timepoint max min from all max/min paths).
% estimate of boundary is fairly accurate assuming animals explored
% entire apparatus and apparatus did not move severely during a time
% point.
df.x_max(1) = 1;
df.y_max(1) = 1;
df.x_min(1) = 1;
df.y_min(1) = 1;
times = unique([df.time_point{:}]);
for idx = times
    temp_idx = [df.time_point{:,1}] == idx;
    df.x_max(temp_idx) = max(x_max(temp_idx,:));
    df.y_max(temp_idx) = max(y_max(temp_idx,:));
    df.x_min(temp_idx) = min(x_min(temp_idx,:));
    df.y_min(temp_idx) = min(y_min(temp_idx,:));
end


end

function df = transform_all_coords(df,dia)
% transforms xy coordinates in dataframe (raw_x and raw_y) that are in
% pixel space into CM. dia is diameter of circular maze.
% Calls FixPos to address potential tracker errors and smooths paths
%
df.x_cm{1} = [];
df.y_cm{1} = [];
% loop through coordinates saved to data table
for i = 1:length(df.file_loc)
    
    x = df.raw_x{i};
    y = df.raw_y{i};
    [ x_cm,y_cm] = transformCoordinates(dia,df.x_max(i),df.x_min(i),df.y_max(i),df.y_min(i),x,y);
    [df.x_cm{i},df.y_cm{i}] = FixPos(x_cm,y_cm,df.ts{i});
    
end

end

function df = locomotion_metrics(df,fr,diameter,binsize)
% computes overall locomotion metrics for all paths in dataframe

% Initalize path variables
df.gen_loco_thigmotaxis(1) = 1;
df.gen_loco_velocity(1) = 1;
df.gen_loco_path_length(1) = 1;
df.gen_loco_search_area(1) = 1;
df.velocity{1} = [];

% loop through coordinates saved to data table
for i = 1:length(df.file_loc)
    
    x = df.x_cm{i};
    y = df.y_cm{i};
    
    % Linear velocity, acceleration, and distance vector
    [velocity, ~, distance_vector] = linear_motion(x,y,fr,2);
    
    % Get index of when the animal is moving greater than 3cm/second
    move = velocity > 3;
    
    % Thigmotaxis
    df.gen_loco_thigmotaxis(i,1) = OF.thigmotaxis(x,y,fr,diameter);
    
    % Search Area
    [~,map] = OF.occ_map(x,y,diameter,binsize,fr);
    df.gen_loco_search_area(i,1) = OF.search_area(map);
    
    % Total distance while moving
    df.gen_loco_path_length(i,1) = sum(distance_vector(move));
    
    % save median velocity while moving
    df.gen_loco_velocity(i,1) = nanmedian(velocity(move));
    
    df.velocity{i,1} = [velocity(1,1);velocity];
    
end

end

function df = stop_metrics(df,fr)
% computes overall locomotion metrics for all paths in dataframe

% Initalize metrics variables
df.stops_num_stops(1)  = nan;
df.stops_median_duration(1)  = nan;
df.stops_median_inter_stop_interval(1)  = nan;
df.timeStopped{1} = nan;
df.stop{1} = nan;
df.inter_stop_interval{1} = nan;
df.stop_start{1} = nan;

% loop through coordinates saved to data table
for i = 1:length(df.file_loc)
    
    x = df.x_cm{i};
    y = df.y_cm{i};
    ts = df.ts{i};
    
    % Linear velocity, acceleration, and distance vecotr
    [velocity, ~, ~] = linear_motion(x,y,fr,0.25);
    
    vel_ts = 0+(1/fr)/2:1/fr:max(ts)-(1/fr)/2;
    
    idx = isnan(velocity);
    velocity(idx) = [];
    vel_ts(idx) = [];
    
    idx = isnan(x) | isnan(y);
    x(idx) = [];
    y(idx) = [];
    ts(idx) = [];
    
    % Find events when rat is stopped (i.e. flip velocity to find 'peaks')
    [events,duration,~] = detect_velocity_epochs(-1*velocity,vel_ts,fr);
    if isnan(events)
        % Save some stop measures for easy plotting later
        df.timeStopped{i} = nan;
        df.stop{i} = nan;
        df.inter_stop_interval{i} = nan;
        df.stop_start{i} = nan;
        
        
        % Calulate some gross run metrics
        df.stops_num_stops(i)  = nan;
        df.stops_median_duration(i)  = nan;
        df.stops_median_inter_stop_interval(i)  = nan;
        continue
    end
    % Grab xy for each stop
    for ii = 1:length(events)
        idx = ts >= events(ii,1) & ts <= events(ii,2);
        stop{ii} = [x(idx),y(idx)];
        stop_start{ii} = events(ii,1);
    end
    
    % Save some stop measures for easy plotting later
    df.timeStopped{i} = duration;
    df.stop{i} = stop;
    df.inter_stop_interval{i} = [NaN;diff(events(:,3))];
    df.stop_start{i} = stop_start;
    
    
    % Calulate some gross run metrics
    df.stops_num_stops(i)  = size(duration,1);
    df.stops_median_duration(i)  = nanmedian(duration);
    df.stops_median_inter_stop_interval(i)  = median(diff(events(:,3)));
    
    clear stop stop_cue_proximity stop_start
end

end

function [events,duration,peak_vel] = detect_velocity_epochs(velocity,vel_ts,fr)
% RH
velocity_thres_sd = 0; % standard deviations
highThresholdFactor = 0.25;
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
    events = nan;
    duration = nan;
    peak_vel = nan;
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
    events = nan;
    duration = nan;
    peak_vel = nan;
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
% figure;
% plot(vel_ts,velocity)
% ys = ylim;
% hold on
% for ep = 1:length(events)
%     %         plot([events(ep,1),events(ep,1)],ylim,'g','LineWidth',3)
%     %         plot([events(ep,2),events(ep,2)],ylim,'r','LineWidth',3)
%     fill([events(ep,1),events(ep,1),events(ep,2),events(ep,2)], [ys(1),ys(2),ys(2),ys(1)], 'm')
% end
% plot(vel_ts,velocity,'k')

%     figure
%     histogram(diff(events(:,3)),200)

end

% Home base distance across days
function primary_hb_dist_mean = hb_dist_by_day(df)
% computes the distance between two homebases based on their centroid.
%
% returns vector of differences same size as length of params. 

% loop through rats 
ratID = unique(df.rat);

% initialize output 
primary_hb_dist_mean = nan(length(ratID),3);

for rat_idx = 1:length(ratID)
    % loop through time point 
    time_idx = 1;
    temp = [df.time_point{:}]; 
    for time_point = unique(temp(df.rat == ratID(rat_idx)))
        % gather coordinates of hb
        hb_day1 = df.HB_center_keep{(df.rat == ratID(rat_idx)) & (df.day == 1) & ([df.time_point{:}]' == time_point)};
        hb_day2 = df.HB_center_keep{(df.rat == ratID(rat_idx)) & (df.day == 2) & ([df.time_point{:}]' == time_point)};
        hb_day3 = df.HB_center_keep{(df.rat == ratID(rat_idx)) & (df.day == 3) & ([df.time_point{:}]' == time_point)};
        
        % check for nan (there was no hb)
        if isnan(hb_day1) | isnan(hb_day2) | isnan(hb_day3)
            primary_hb_dist_mean(rat_idx,time_idx)  = nan;
            continue
        
        end
        
        % compute distance vector between coordinates 
        primary_hb_dist_1v2 = sqrt((hb_day1(1,1) - hb_day2(1,1))^2 + (hb_day1(1,2) - hb_day2(1,2))^2);
        primary_hb_dist_1v3 = sqrt((hb_day1(1,1) - hb_day3(1,1))^2 + (hb_day1(1,2) - hb_day3(1,2))^2);
        primary_hb_dist_2v3 = sqrt((hb_day2(1,1) - hb_day3(1,1))^2 + (hb_day2(1,2) - hb_day3(1,2))^2);
        
        primary_hb_dist_mean(rat_idx,time_idx) = nanmean([primary_hb_dist_1v2,primary_hb_dist_1v3,primary_hb_dist_2v3]);
        
        time_idx = time_idx+1;
        
        clear primary_hb_dist_1v2 primary_hb_dist_1v3 primary_hb_dist_2v3
    end
    
end
% reorganize to long format 
primary_hb_dist_mean = [ratID repmat(4,length(ratID),1) primary_hb_dist_mean(:,1); ...
    ratID repmat(7,length(ratID),1) primary_hb_dist_mean(:,2); ...
    ratID repmat(12,length(ratID),1) primary_hb_dist_mean(:,3);];
end

function df = home_base_analysis(df,fr,binsize,diameter,upFactor)
% identifies homebase using kmeans of high occupancy coordinates and
% obtains metrics relating to stopping behaviors.

% initalize variables to save in dataframe
% Total duraiton
df.hb_duration(1) = nan;

% Distance between homebase and stops
df.hb_stop_dist(1) = nan;

% Number of stops that are within 25cm of home base
df.hb_close_stop(1) = nan;

% Total number of HB stops
df.hb_home_base_stops(1) = nan;

% First time to home base
df.hb_time2HB_keep(1) = nan;

% number of home bases
df.hb_HB_count_keep(1) = nan;

% testing home base area 
df.hb_area{1} = [];

% Total duration all 
df.hb_duration_all{1} = [];

% Keep the home base center for between day tests
df.HB_center_keep{1} = [];

% hb min distance from wall
df.hb_min_dist_from_wall(1) = nan;

for i = 1:length(df.file_loc)
    disp(num2str(i))
    
    if isempty(df.x_cm{i})
        continue
    end
    
    x = df.x_cm{i};
    y = df.y_cm{i};
    ts = df.ts{i};
    
    % Linear velocity, acceleration, and distance vecotr
    [velocity, ~, ~] = linear_motion(x,y,fr,0.25);
    vel_ts = 0+(1/fr)/2:1/fr:max(ts)-(1/fr)/2;
    
    velocity = [velocity(1); velocity];
    
    idx = isnan(velocity);
    velocity(idx) = [];
    vel_ts(idx) = [];
    x(idx) = [];
    y(idx) = [];
    ts(idx) = [];
    
    % Find events when rat is stopped (i.e. flip velocity to find 'peaks')
    [events,~,~] = detect_velocity_epochs(-1*velocity,vel_ts,fr);
    x_stop = [];
    y_stop = [];
    ts_stop = [];
    
    % if there are no moving slow events, then there are no homebases.
    if isnan(events)
            % Total duraiton
            df.hb_duration(i) = nan;
            
            % Total duraiton all hb 
            df.hb_duration_all{i} = [];
            
            % Distance between homebase and stops
            df.hb_stop_dist(i) = nan;
            
            % Number of stops that are within 25cm of home base
            df.hb_close_stop(i) = nan;
            
            % Total number of HB stops
            df.hb_home_base_stops(i) = nan;
            
            % First time to home base
            df.hb_time2HB_keep(i) = nan;
            
            % number of home bases
            df.hb_HB_count_keep(i) = nan;
            
            % hb center
            df.HB_center_keep{i} = nan;
            
            % hb min distance from wall
            df.hb_min_dist_from_wall(i) = nan;
            
            continue
    end
    
    for e = 1:length(events)
        idx = ts >= events(e,1) & ts <= events(e,2);
        try
        [ ~,~, stop_center{e}] = min_encl_ellipsoid(x(idx),y(idx));
        catch
            stop_center{e} = [0;0];
        end
        x_stop = [x_stop; x(idx)];
        y_stop = [y_stop; y(idx)];
        ts_stop = [ts_stop; ts(idx)];
    end
    
    % Occupancy Map Measures (search area, home bases, home base area, all home base measures)
    [~,map] = OF.occ_map(x,y,diameter,binsize,fr);
    
    % High occupancy coordinates
    [~,~,home_base_x,home_base_y,df.hb_area{i},upHBmap] = segmentImage('map',map,'upscalefac',upFactor,'figs',false); %Store area
    
    %% Preprocess high occupancy boundaries
    
    for hb = 1:length(home_base_x)
        
        % rescale coordinates back to pool size
        x_home = rescale([home_base_x{hb},upFactor+1,size(upHBmap,2)-upFactor],-(diameter/2),(diameter/2));
        y_home = rescale([home_base_y{hb},upFactor+1,size(upHBmap,2)-upFactor],-(diameter/2),(diameter/2));
        
        % logical index for all coordinates inside home base
        in_home = inpolygon(x,y,x_home(1:end-2)',y_home(1:end-2)');
        
        % index for frames inside home base for at least 2seconds
        out_home = contiguousframes(in_home,fr*2); %has to be inside of hb for at least 2 sec to count as entry
        
        % Find proportion of time moving slow in home base
        slow_in_homebase{hb} = sum(velocity(out_home(1:end-1)) < 3)/(nansum(out_home));
        
        HB_duration{hb} = sum(out_home)/fr;
            
        if slow_in_homebase{hb} < .70

            HB_stop_dist{hb} = NaN;
            HB_close_stop{hb} = NaN;
            home_base_stops{hb} = NaN;
            time2HB{hb} = NaN;
            HB_duration{hb} = NaN;
            HB_center_keep{hb} = nan;
            hb_wall_proximity{hb} = nan;
            
        else
            % boundary of home base
            HBBound = [x_home(1,1:end-2)',y_home(1,1:end-2)'];
            
            % Calculate the center of the home base
            try
               [ ~,~, tempC] = min_encl_ellipsoid(HBBound(:,1),HBBound(:,2));
               HBcenter= [tempC(1,1),tempC(2,1)];

            catch
                % if it fails, just take first hb coordinate 
                HBcenter= [HBBound(1,1),HBBound(1,2)];
            end

            [HB_stop_dist{hb},HB_close_stop{hb}] = stops_to_homebase(HBcenter,stop_center);
            
            [home_base_stops{hb},time2HB{hb}] = home_base_metics(x_home,y_home,x_stop,y_stop,stop_center,ts_stop);
            
            HB_center_keep{hb} = HBcenter;
            
            hb_wall_proximity{hb} = hb_proximity_to_wall(x_home,y_home,diameter);
        end
        
    end
    % save overall metrics and primary home base metrics in df
    
    % Unpack HB coord and keep max given moving slow
    slow_in_hb_idx = [slow_in_homebase{:}] > .70;
    
    [~,duration_idx] = max([HB_duration{:}]);
    
    % Total duraiton
    df.hb_duration(i) = HB_duration{duration_idx};
    df.hb_duration_all{i} = [HB_duration{:}];
    
    % keep center of primary home base
    df.HB_center_keep{i} = HB_center_keep{duration_idx};
    
    % hb min distance from wall
    df.hb_min_dist_from_wall(i) = hb_wall_proximity{duration_idx};

    % Distance between homebase and stops
    df.hb_stop_dist(i) = HB_stop_dist{duration_idx};
    
    % Number of stops that are within 25cm of home base
    df.hb_close_stop(i) = HB_close_stop{duration_idx} ;
    
    % Total number of HB stops
    df.hb_home_base_stops(i) = home_base_stops{duration_idx} ;
    
    % First time to home base
    df.hb_time2HB_keep(i) = time2HB{duration_idx};
    
    % number of home bases
    df.hb_HB_count_keep(i) = sum(slow_in_hb_idx);
    
    clear HB_duration HB_stop_dist HB_close_stop home_base_stops time2HB slow_in_homebase stop_center
    
end


end
function [HB_stop_dist,HB_close_stop] = stops_to_homebase(HBcenter,stop_center)

% Average Proximity of stops from hb center
disttemp = zeros(size(stop_center,2),1);
for i = 1:size(stop_center,2)
    temp = stop_center{i};
    dist = sqrt((HBcenter(1,1)-temp(1,1)).^2+(HBcenter(1,2)-temp(2,1)).^2);
    disttemp(i,1) = nanmean(dist);
end

HB_stop_dist = nanmean(disttemp); %average stop distance
HB_close_stop = sum(disttemp<25); %number of stops within 25cm from hb center
% HB_stop_dist_vector = disttemp;

end

function hb_wall_proximity = hb_proximity_to_wall(x_home,y_home,diameter)

% boundary of maze (assumes circular maze)
boundary_coords = [[sin(0:pi/360:2*pi)*(diameter/2)]',[cos(0:pi/360:2*pi)*(diameter/2)]'];

% distance from cue boundary
for r = 1:length(x_home)
    distances = sqrt(sum(bsxfun(@minus, boundary_coords, [x_home(r),y_home(r)]).^2,2));
    move_cue_distance(r) = unique(distances(distances==min(distances))); %Find minimum distance from cue boundary to hb center
end

% find minimum distance
hb_wall_proximity = nanmean(move_cue_distance);

end

function [home_base_stops,time2HB] = home_base_metics(x_home,y_home,x_stop,y_stop,stop_center,ts_stop)

% This finds stops that occur in the home base boundary
temp_center = horzcat(stop_center{:})';
home_base_stops = sum(inpolygon(temp_center(:,1),temp_center(:,2),x_home(1:end-2)',y_home(1:end-2)'));

% index for the time it takes to reach the home base
time2HB_idx = inpolygon(x_stop,y_stop,x_home(1:end-2)',y_home(1:end-2)');
time2HB = ts_stop(time2HB_idx);
if isempty(time2HB)
    time2HB = nan;
else
    time2HB = time2HB(1);
end

% shortest distance of home base to wall


end