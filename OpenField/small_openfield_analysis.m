% small_openfield_analysis examines locomotion and home base behaviors of
% TgF344-AD and F344 controls in small open field (76.5cm) in 10min
% sessions across days (days 1,2,3) and time-points (4 months, 7 months, 10
% months). Positional data was tracked using anymaze and xy coordinates
% were saved as .csv files using the Process_SLK.ipynb 

% initalize data table
df = table;

% Import CSV files 
files = dir("H:/OF/coords_csv/**/*.csv");

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
%% Preprocess paths (smooth) 
for i = 1:length(df.file_loc)
    temp = readtable(df.file_loc{i});
    x = temp.x;
    y = temp.y;
    % transform coordinates to cm 
    [params.noseCM{i}(:,1),params.noseCM{i}(:,2)] = transformCoordinates...
        (76.5,params.Xmax{i},params.Xmin{i},params.Ymax{i},params.Ymin{i},...
        x,y); %Converts the coordinates into CM
end

%% Locomotion (path length, running speed, occupancy)

%% Stops (stop frequency, mean/median stop duration) 

%% Homebase (number of home bases) 

%% Save outcome measures to csv for analysis in R 

