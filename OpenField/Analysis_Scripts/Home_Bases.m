%% Home Bases
% Examines high occupancy places in the environment

% Stops - epochs of time where the animal moves <3cm/s for at least 2
% seconds

% General maze exploration
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/'; % maxmin.mat, cueCoords.mat

% load data file
load('/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/params.mat'); % maxmin.mat, cueCoords.mat

% Set some measurement parameters
fr = 30; % 30 hz sample rate
binsize = 3;
diameter = params.dia{1};
upFactor = 15; %used in segmentimage - this upscales the heatmap to get a smooth boundary


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
    
    velocity = [velocity(1); velocity];
    
    idx = isnan(velocity);
    velocity(idx) = [];
    vel_ts(idx) = [];
    x(idx) = [];
    y(idx) = [];
    ts(idx) = [];
    
    % Find events when rat is stopped (i.e. flip velocity to find 'peaks')
    [events,duration,~] = detect_velocity_epochs(-1*velocity,vel_ts,fr);
    x_stop = [];
    y_stop = [];
    ts_stop = [];
    for e = 1:length(events)
        idx = ts >= events(e,1) & ts <= events(e,2);
        [ ~,~, stop_center{e}] = min_encl_ellipsoid(x(idx),y(idx));
        x_stop = [x_stop; x(idx)];
        y_stop = [y_stop; y(idx)];
        ts_stop = [ts_stop; ts(idx)];
    end
    
    % Occupancy Map Measures (search area, home bases, home base area, all home base measures)
    [~,map] = OF.occ_map(x,y,diameter,binsize,fr);
    
    % High occupancy coordinates
    [~,~,home_base_x,home_base_y,fieldarea,upHBmap] = segmentImage('map',map,'upscalefac',upFactor); %Store area
    
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
            HB_avg_dist{hb} = NaN;
            HB_max_dist{hb} = NaN;
            HB_min_dist{hb} = NaN;
            HBdist2Cue{hb} = NaN;
            HB_stop_dist{hb} = NaN;
            HB_close_stop{hb} = NaN;
            home_base_stops{hb} = NaN;
            time2HB{hb} = NaN;
            HB_duration{hb} = NaN;
            
        else
            % boundary of home base
            HBBound = [x_home(1:end-2)',y_home(1:end-2)'];
            
            % Calculate the center of the home base
            [ ~,~, tempC] = min_encl_ellipsoid(HBBound(:,1),HBBound(:,2));
            HBcenter= [tempC(1,1),tempC(2,1)];
            
            [HB_avg_dist{hb},~,~] = distance_between_homebase(HBcenter);
            
            if isnan(cue_cords)
                HBdist2Cue{hb} =  homebase_dist_from_cue(params.cueCM{i - 1},HBcenter);
            else
                HBdist2Cue{hb} =  homebase_dist_from_cue(cue_cords,HBcenter);
                
            end
            
            [HB_stop_dist{hb},HB_close_stop{hb}] = stops_to_homebase(HBcenter,stop_center);
            
            [home_base_stops{hb},time2HB{hb}] = home_base_metics(x_home,y_home,x_stop,y_stop,stop_center,ts_stop);
            
        end
        
    end
    
    % Unpack HB coord and keep max given moving slow
    slow_in_hb_idx = [slow_in_homebase{:}] > .70;
   
    [~,duration_idx] = max([HB_duration{:}]);
   
    % Total duraiton
    HB_duration_keep(i) = HB_duration{duration_idx};
    
    % Distance between other home bases
    HB_avg_dist_keep(i) = HB_avg_dist{duration_idx};
    
    % Distance between homebase and cue
    HBdist2Cue_keep(i) = HBdist2Cue{duration_idx};
    
    % Distance between homebase and stops 
    HB_stop_dist_keep(i) = HB_stop_dist{duration_idx};
    
    % Number of stops that are within 25cm of home base
    HB_close_stop_keep(i) = HB_close_stop{duration_idx} ;
    
    % Total number of HB stops 
    home_base_stops_keep(i) = home_base_stops{duration_idx} ;
    
    % First time to home base
    time2HB_keep(i) = time2HB{duration_idx};
    
    % number of home bases 
    HB_count_keep(i) = sum(slow_in_hb_idx);

    clear slow_in_homebase time2HB home_base_stops ...
        HB_close_stop HB_stop_dist HBdist2Cue HB_min_dist HB_max_dist...
        HB_avg_dist HB_duration stop_center home_base_x home_base_y

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
vars = {'subID','group','day',...
    'duration',...
    'avg_dist_between_HB',...
    'distance_to_cue',...
    'stop_distance',...
    'num_close_stops',...
    'num_stops',...
    'time_to_HB',...
    'num_of_HB'...
    };

df = table(subID,group,day,HB_duration_keep',HB_avg_dist_keep',...
    HBdist2Cue_keep',HB_stop_dist_keep',...
    HB_close_stop_keep',home_base_stops_keep',time2HB_keep',HB_count_keep','VariableNames',vars);


% Writetable
writetable(df,[save_path,'home_bases.csv'])


%% Local functions

% Movement proximity to cue
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

% Discard ripples that are too short
events(duration < min_event_duration,:) = NaN;
events = events((all((~isnan(events)),2)),:);

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

% Home base distance from cue
function HBdist2Cue =  homebase_dist_from_cue(cueCM,HBcenter)
% calculating proximity of high occupancy coordinates center from the
% cue boundary

k = convhull(cueCM(:,1),cueCM(:,2));
cueBoundary = [cueCM(k,1),cueCM(k,2)];


distances = sqrt(sum(bsxfun(@minus, cueBoundary, [HBcenter(1,1),HBcenter(1,2)]).^2,2));
HBdist2Cue = unique(distances(distances == min(distances))); %Find minimum distance from cue boundary to hb center


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


function [HB_avg_dist,HB_max_dist,HB_min_dist] = distance_between_homebase(HBcenter)
% input:
%   - HBcenter: vector of HB centers
% output:
%   -HB_avg_dist: average distance between home bases
%   -HB_max_dist: maximum distance between home bases
%   -HB_min_dist: minimum distance between home bases
%
% Calculate distance measures between high occupancy coordinates centers
temp = [];
if size(HBcenter,2) > 1
    for hb = 1:size(HBcenter,2)
        temp = [temp; HBcenter(1,hb)];
    end
    HB_avg_dist = nanmean(pdist(temp));
    HB_max_dist = max(pdist(temp));
    HB_min_dist = min(pdist(temp));
else
    HB_avg_dist = NaN;
    HB_max_dist = NaN;
    HB_min_dist = NaN;
end
end

function [home_base_stops,time2HB] = home_base_metics(x_home,y_home,x_stop,y_stop,stop_center,ts_stop)

% This finds stops that occur in the home base boundary
temp_center = horzcat(stop_center{:})';
home_base_stops = sum(inpolygon(temp_center(:,1),temp_center(:,2),x_home(1:end-2)',y_home(1:end-2)'));

% index for the time it takes to reach the home base
time2HB_idx = inpolygon(x_stop,y_stop,x_home(1:end-2)',y_home(1:end-2)');
time2HB = ts_stop(time2HB_idx);
time2HB = time2HB(1);

end