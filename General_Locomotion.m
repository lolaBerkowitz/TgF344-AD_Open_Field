
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
    
    % Linear velocity, acceleration, and distance vecotr
    [velocity, acceleration, distance_vector] = linear_motion(x,y,fr,3);
    
    % Get index of when the animal is moving greater than 3cm/second
    move = velocity > 3;
    
    % Thigmotaxis 
    out(i,1) = OF.thigmotaxis(x,y,fr,diameter);

    
    % Search Area
    [c,map] = OF.occ_map(x,y,diameter,binsize,fr);
    search_area(i,1) = OF.search_area(map);
    
    % Head Direction (use head and neck to compute head direction)
    theta = wrapTo360(rad2deg(atan2(params.neckCM{1,1}(:,2)-params.headCM{1,1}(:,2),...
        params.neckCM{1,1}(:,1)-params.headCM{1,1}(:,1))));
    
    % Angular velocty
    angular_velocity = insta_angvel(deg2rad(theta),fr,3);
    
    % Compile median values
    median_velocity(i,1) = nanmedian(velocity(move));
    median_acceleration(i,1) = nanmedian(abs(acceleration(move)));
    median_angular_velocity(i,1) = nanmedian(abs(angular_velocity(move)));
    
    % Total distance
    path_length(i,1) = sum(distance_vector(move));
    
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
    'search_area','thigmotaxis'};

df = table(subID,group,day,path_length,median_velocity,median_acceleration,...
    median_angular_velocity,search_area,out,'VariableNames',vars);


% Writetable 
writetable(df,[save_path,'general_locomotion.csv'])

