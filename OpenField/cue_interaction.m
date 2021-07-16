%% Interaction with cue 
% run after movement_epochs and stop_epochs 

% General maze exploration
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/'; % maxmin.mat, cueCoords.mat

% load data file
load('/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/params.mat'); % maxmin.mat, cueCoords.mat

params.time_2_cue{1} = [];
%% Loop through subjects
for i = 1:length(params.subID)
    
    % Compute measures from back coordinates
    x = params.backCM{i}(:,1);
    y = params.backCM{i}(:,2);
    head_x = params.headCM{i}(:,1);
    head_y = params.headCM{i}(:,2);
    neck_x = params.neckCM{i}(:,1);
    neck_y = params.neckCM{i}(:,2);
    ts = params.ts{i}(:,1);
    cue_cords = params.cueCM{i};
    
    
    %Compute circuity and cue proximity, between run segments
    if isnan(cue_cords)
        % Check when neck is in cue area
        params.time_2_cue{i} = time_to_cue(neck_x,neck_y,ts,params.cueCM{i - 1});
        time_2_cue(i) = params.time_2_cue{i};
    else
        % Check when neck is in cue area
        params.time_2_cue{i} = time_to_cue(neck_x,neck_y,ts,cue_cords);
        time_2_cue(i) = params.time_2_cue{i};
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

df = table(subID,group,day,time_2_cue','VariableNames',{'subID','group','day','time_2_cue'});


% Writetable
writetable(df,[save_path,'time_2_cue.csv'])




function first_visit_time = time_to_cue(x,y,ts,cue_cords)
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

% Distances within 5cm of cue boundary 

idx = move_cue_distance < 10;
if sum(idx) == 0
    first_visit_time = NaN;
else
    in_cue_area = ts(idx);
    first_visit_time = in_cue_area(1);
end

end
