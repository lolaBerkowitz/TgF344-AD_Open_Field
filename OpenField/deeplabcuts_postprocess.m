% deeplabcuts_postprocess RH May 2019, edits by LB June 2019

figs = 0;

% Hardcoded parameters 
dia = 202 ; %diameter of maze in cm

%Initialize path to DLC output
path_to_files ='/Users/lauraberkowitz/Google Drive/Manuscripts/In Progress/TgF344-AD_OF/Data/Tracked'; %contains DLC csv
path_to_data = '/Users/lauraberkowitz/Google Drive/Manuscripts/In Progress/TgF344-AD_OF/Data/'; % maxmin.mat, cueCoords.mat

files = struct2table(dir([path_to_files,filesep,'**/*.csv']));

%Initialize data table and variables
params=table;
params.subID{1}=[];
params.head{1}=[];
params.neck{1}=[];
params.back{1}=[];
params.tail{1}=[];

%% Compile environment Max/Min for transformation into cm and cue Coords for
%computing cue related measures in OF_postprocess

load([path_to_data,'maxmin.mat']);

for i = 1:length(files.name)
    
    params.subID{i} = extractBefore(files.name{i},'DLC');
    params.dia{i} = dia;
    % load header
    fileID = fopen(fullfile(files.folder{i},files.name{i}),'r');
    dataArray = textscan(fileID, '%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]',...
        3-2+1, 'Delimiter',',', 'TextType', 'string', 'HeaderLines',...
        2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    header = [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
    
    % load data
    fileID = fopen(fullfile(files.folder{i},files.name{i}),'r');
    dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]',...
        'Delimiter', ',','TextType', 'string', 'EmptyValue', NaN,...
        'HeaderLines' ,4-1,'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    tsxy= [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
    
    % load coded behaviors(obtained from Behavior_code.m)
    if exist([files.folder{i},filesep,params.subID{i},'.mat'])
        behav = load([files.folder{i},filesep,params.subID{i},'.mat']);
    else
        behav = nan;
    end
    
    % convert ts into seconds
    tsxy(:,1)=0:1/30:tsxy(end,1)/30;
    
%     Create outer limit to remove points outside the maze. 
    xunit=sin(0:pi/360:2*pi)*((maxmin.xmax{i}-maxmin.xmin{i})/2)*1.1+median([maxmin.xmax{i},maxmin.xmin{i}]);
    yunit=cos(0:pi/360:2*pi)*((maxmin.ymax{i}-maxmin.ymin{i})/2)*1.1+median([maxmin.ymax{i},maxmin.ymin{i}]);

    % Turn low liklihood point into NAN and smooth remaining coords
    idx=tsxy(:,contains(header(2,:),'likelihood')) < .80;
    xloc=find(contains(header(2,:),'x'));
    yloc=find(contains(header(2,:),'y'));
    tstart=behav.Trial_start(end,1);
    
    tend = tstart + (30*60); %30=30min trial duration, 60=seconds/min, 
    for l=1:size(idx,2)
        tsxy(idx(:,l),xloc(l):yloc(l)) = NaN;
        
        tsxy(~inpolygon(tsxy(:,xloc(l)),tsxy(:,yloc(l)),xunit,yunit),xloc(l):yloc(l))=NaN;
        
        [tsxy(:,xloc(l)),tsxy(:,yloc(l))] = FixPos(tsxy(:,xloc(l)),tsxy(:,yloc(l)),tsxy(:,1));
    end
    
    %Keep only points that are within the trial start time plus 30min
    tsxy=tsxy(tsxy(:,1)>=tstart & tsxy(:,1)<=tend,:);
    
    %Reset timestamps to zero 
    tsxy(:,1)=(tsxy(:,1)-tsxy(1,1));
    
    % Save data
    params.ts{i}=tsxy(:,1);
    params.head{i}(:,1) = tsxy(:,2);
    params.head{i}(:,2) = tsxy(:,3);
    params.neck{i}(:,1)=tsxy(:,5);
    params.neck{i}(:,2)=tsxy(:,6);
    params.back{i}(:,1) = tsxy(:,8);
    params.back{i}(:,2) = tsxy(:,9);
    params.tail{i}(:,1)=tsxy(:,11);
    params.tail{i}(:,2)=tsxy(:,12);
    
    if figs==1
        figure;
        plot(params.nose{i}(:,1),params.nose{i}(:,2),'.k')
        hold on
        nPlot = plot(NaN,NaN,'ro','MarkerFaceColor','r');
        hPlot = plot(NaN,NaN,'bo','MarkerFaceColor','b');
        bPlot = plot(NaN,NaN,'go','MarkerFaceColor','g');
        buttPlot = plot(NaN,NaN,'co','MarkerFaceColor','c');
        
        
        xlim([min(params.nose{i}(:,1)) max(params.nose{i}(:,1))]);
        ylim([min(params.nose{i}(:,2)) max(params.nose{i}(:,2))]);
        % iterate through each point on line
        for k=1:length(params.nose{i}(:,1))
            set(nPlot,'XData',params.head{i}(k,1),'YData',params.head{i}(k,2));
            set(hPlot,'XData',params.neck{i}(k,1),'YData',params.neck{i}(k,2));
            set(bPlot,'XData',params.back{i}(k,1),'YData',params.back{i}(k,2));
            set(buttPlot,'XData',params.tail{i}(k,1),'YData',params.tail{i}(k,2));
            
            pause(1/240);
        end
        close
    end
end

load([path_to_data,filesep,'cueCoords.mat']);

params.Xmax = maxmin.xmax;
params.Xmin = maxmin.xmin;
params.Ymax = maxmin.ymax;
params.Ymin = maxmin.ymin;
params.cueCoords = cueCoords.coords;

save([path_to_data,'params_raw_data.mat'])
