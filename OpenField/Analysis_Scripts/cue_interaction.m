%% Interaction with cue 
% run after movement_epochs and stop_epochs 

% LB 08/2021

% General maze exploration
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/'; % maxmin.mat, cueCoords.mat

% load data file
load('/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/params.mat'); % maxmin.mat, cueCoords.mat

% video parameters 
fr = 30; % frame rate of video is 30Hz

params.time_2_cue{1} = [];
params.cue_close_stops{1} = [];
params.prop_time_in_cue_area{1}=[];
%% Loop through subjects
for i = 1:length(params.subID)
    
    % Compute measures from back coordinates
    x = params.backCM{i}(:,1);
    y = params.backCM{i}(:,2);
    ts = params.ts{i}(:,1);
    cue_cords = params.cueCM{i};
    stops = params.stop{i};
    
    
    %Compute cue interaction metrics (time to cue, stops near cue, time
    % spent in cue area)
    if isnan(cue_cords)
        % Check when back is in cue area
        [temp_time_in_cue_area,first_visit_time] = time_to_cue(x,y,ts,params.cueCM{i - 1},fr); 
        params.time_2_cue{i} = first_visit_time;
        params.prop_time_in_cue_area{i} = temp_time_in_cue_area;
        time_2_cue(i) = first_visit_time;
        time_in_cue_area(i) = temp_time_in_cue_area;
        
        % compute number of close stops 
        [cue_close_stop,cue_stop_dist] = stops_in_cue_zone(stops,params.cueCM{i - 1});
        params.cue_close_stops{i} = cue_stop_dist;
        close_stop(i) = cue_close_stop;
    else
        % Check when back is in cue area
        [temp_time_in_cue_area,first_visit_time] = time_to_cue(x,y,ts,cue_cords,fr);
        params.time_2_cue{i} = first_visit_time;
        params.prop_time_in_cue_area{i} = temp_time_in_cue_area;
        time_2_cue(i) = first_visit_time;
        time_in_cue_area(i) = temp_time_in_cue_area;

        
        
         % compute number of close stops 
        [cue_close_stop,cue_stop_dist] = stops_in_cue_zone(stops,cue_cords);
        params.cue_close_stops{i} = cue_stop_dist;
        close_stop(i) = cue_close_stop;

    end
    
    

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

%% Save to google drive

df = table(subID,group,day,time_2_cue',time_in_cue_area','VariableNames',{'subID','group','day','time_2_cue','time_in_cue_area'});


% Writetable
writetable(df,[save_path,'time_2_cue.csv'])


function [cue_close_stop,cue_stop_dist] = stops_in_cue_zone(stops,cue_cords)
% for a given set of stops (cell array of stops from params) and cue coordinates (boundary coordinates of cue)
%, this function will return the number of stops that were close to cue and the distances for each
% stop to cue boudary. 
%

cue_x = cue_cords(:,1);
cue_y = cue_cords(:,2);

% convex hull to close shape of cue
k = convhull(cue_x,cue_y);
cue_boundary = [cue_x(k) cue_y(k)];

% Average Proximity of stops from hb center
stop_center = zeros(size(stops,2),1);
for i = 1:size(stop_center,1)
    stop_x = stops{1,i}(:,1);
    stop_y = stops{1,i}(:,2);
    for ii = 1:length(stop_x)
            % distance from each point in stop to cue boundary
            disttemp = sqrt(sum(bsxfun(@minus, cue_boundary, [stop_x(ii),stop_y(ii)]).^2,2));
    end
    
    % compile minimum distance of each stop to cue boundary
    cue_stop_dist(i,1) = nanmin(disttemp);
end

% count number of stops that were within 25cm from cue boundary
cue_close_stop = sum(cue_stop_dist<25); 
end

function [time_in_cue_area,first_visit_time] = time_to_cue(x,y,ts,cue_cords,fr)
cue_x = cue_cords(:,1);
cue_y = cue_cords(:,2);

% convex hull to close shape of cue
k = convhull(cue_x,cue_y);
cue_boundary = [cue_x(k) cue_y(k)];

% distance from cue boundary
for r = 1:length(x)
    distances = sqrt(sum(bsxfun(@minus, cue_boundary, [x(r),y(r)]).^2,2));
    if isnan(x(r)) || isnan(y(r))
        move_cue_distance(r) = NaN;
    else
        move_cue_distance(r) = unique(distances(distances==min(distances))); %Find minimum distance from cue boundary to coordinate

    end
    
end

% Distances within 25cm of cue boundary 

idx = move_cue_distance < 25;
if sum(idx) == 0
    time_in_cue_area = 0;
    first_visit_time = NaN;
else
    in_cue_area = ts(idx);
    time_in_cue_area = length(in_cue_area)/fr; % to save in seconds
    first_visit_time = in_cue_area(1);
end

end
