
%%
%%%%%%%%%%% Occupancy Map for Day 1 (Figure 1)
% Fig 1 example paths plotted using plotOFPaths.m
% Fig 1 statistical test plots done using dabest in python
% TgF344_OF_figs_v2.py script

% load('params_V8');
load('/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/params.mat'); % maxmin.mat, cueCoords.mat
param_idx=params.subID;


%% General Locomotion - Example paths for both groups on day 1
% Day 2 all paths per group. Save to NPY for plotting in python. 

params_tg = params(contains(params.subID,'Tg') & contains(params.subID,'D2'),:);

paths_tg = [];
for i = 1:length(params_tg.subID)
    
    paths_tg = [paths_tg;params_tg.backCM{i}]; 
    
end

writeNPY(paths_tg, '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/all_paths_tg_day2.npy')

params_wt = params(contains(params.subID,'WT') & contains(params.subID,'D2'),:);

paths_wt = [];
for i = 1:length(params_wt.subID)
    
    paths_wt = [paths_wt;params_wt.backCM{i}]; 
    
end

writeNPY(paths_wt, '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/all_paths_wt_day2.npy')

% Day 1 all paths per group. Save to NPY for plotting in python. 

params_tg = params(contains(params.subID,'Tg') & contains(params.subID,'D1'),:);

paths_tg = [];
for i = 1:length(params_tg.subID)
    
    paths_tg = [paths_tg;params_tg.backCM{i}]; 
    
end

writeNPY(paths_tg, '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/all_paths_tg.npy')

params_wt = params(contains(params.subID,'WT') & contains(params.subID,'D1'),:);

paths_wt = [];
for i = 1:length(params_wt.subID)
    
    paths_wt = [paths_wt;params_wt.backCM{i}]; 
    
end

writeNPY(paths_wt, '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/all_paths_wt.npy')

%%
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/figs/MatlabFigs/';
% set figure defaults
fig=figure('DefaultAxesFontSize',8,'defaultAxesFontName','Serif','defaultTextFontName','Serif');
fig.Color = [1,1,1];
[fig_width_in, fig_height_in] = set_size('thesis', .75, [1,1]);
set(fig,'Position',[835 270 fig_width_in fig_height_in])
plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k')
hold on; plot((cos(linspace(-pi,pi,1000))+0)*(101*.8),(sin(linspace(-pi,pi,1000))+0)*(101*.8),'--k')
p1 = plot(paths_tg(:,1),paths_tg(:,2),'LineWidth',1,'Color','#601a4a');
p1.Color(4) = 0.5;
axis image
axis off

saveas(fig,[save_path,filesep,'Tg_D1_paths','.pdf'],'pdf')


% set figure defaults
fig=figure('DefaultAxesFontSize',8,'defaultAxesFontName','Serif','defaultTextFontName','Serif');
fig.Color = [1,1,1];
[fig_width_in, fig_height_in] = set_size('thesis', .75, [1,1]);
set(fig,'Position',[835 270 fig_width_in fig_height_in])
plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k')
hold on; plot((cos(linspace(-pi,pi,1000))+0)*(101*.8),(sin(linspace(-pi,pi,1000))+0)*(101*.8),'--k')
p1 = plot(paths_wt(:,1),paths_wt(:,2),'LineWidth',1,'Color','#9c9eb5');
p1.Color(4) = 0.5;
axis image
axis off

saveas(fig,[save_path,filesep,'WT_D1_paths','.pdf'],'pdf')


%% Example Segments for One Tg and one Wt animal (Figure 2)
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/figs/MatlabFigs/path_circuity/';

%create colormap the size of the unique values in data -
% now each line can be colorcoded on the same scale by the circuity value
temp_circ = [];
for all = 1:length(param_idx)
    temp_circ = [temp_circ params.run_circ{all}];
end

unique_circ = unique([temp_circ],'sorted');
cm_magma = magma(size(unique_circ,2));

for rat = 1:length(param_idx)
    fig=figure;
    fig.Color = [1 1 1];
    [fig_width_in, fig_height_in] = set_size('thesis', .75, [1,1]);
    set(fig,'Position',[835 270 fig_width_in fig_height_in])
    plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k'); hold on;
    
    title(extractBefore(param_idx(rat,1),'_'))
    axis off
    %     for i = 1:size(params.runs{rat},2)
    path_circ = params.run_circ{rat};
    [~,path_idx] = min(path_circ);
    temp_run = params.runs{rat}{1,path_idx};
    path_circ = params.run_circ{rat}(path_idx);
    idx = ismember(unique_circ',path_circ');
    plot(temp_run(:,1),temp_run(:,2),'Color',cm_magma(idx,:),'LineWidth',3); hold on;
    %     end
    axis image
    ax = axes;
    c = colorbar(ax);
    ax.Visible = 'off';
    colormap(cm_magma)
    
    saveas(fig,[save_path,filesep,param_idx{rat,1},'_circ_paths_min','.svg'],'svg')
    
    close all
    
end

%%
%Tg examples 3(min),5(max)
% Wt examples 45(min),29(max)

Tg1_runs = params.runs{3};
for i = 1:size(Tg1_runs,2)
    temp_run = Tg1_runs{1,i};
    temp_run(isnan(temp_run(:,1)),:)=[];
    Tg1_circ(i,1)=circuity(temp_run(:,1),temp_run(:,2));
end
Tg2_runs = params.runs{9};
for i = 1:size(Tg2_runs,2)
    temp_run = Tg2_runs{1,i};
    temp_run(isnan(temp_run(:,1)),:)=[];
    Tg2_circ(i,1)=circuity(temp_run(:,1),temp_run(:,2));
end
Wt1_runs = params.runs{45};
for i = 1:size(Wt1_runs,2)
    temp_run = Wt1_runs{1,i};
    temp_run(isnan(temp_run(:,1)),:)=[];
    Wt1_circ(i,1)=circuity(temp_run(:,1),temp_run(:,2));
end
Wt2_runs = params.runs{29};
for i = 1:size(Wt2_runs,2)
    temp_run = Wt2_runs{1,i};
    temp_run(isnan(temp_run(:,1)),:)=[];
    Wt2_circ(i,1)=circuity(temp_run(:,1),temp_run(:,2));
end


% Find unique values to pull corresponding color to colormap - sort to
% maintain value relationships
unique_circ = unique([Tg1_circ;Wt1_circ;Tg2_circ;Wt2_circ],'sorted');

fig=figure;
fig.Color = [1 1 1];
subaxis(2,2,2)
plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k'); hold on;
title(extractBefore(param_idx(3,1),'_'))
axis off
for i = 1:size(Tg1_runs,2)
    temp_run = Tg1_runs{1,i};
    temp_circ = Tg1_circ(i,1);
    ismember(unique_circ,Tg1_circ(1,1));
    plot(temp_run(:,1),temp_run(:,2),'Color',cm_magma(ismember(unique_circ,Tg1_circ(i,1)),:),'LineWidth',3); hold on;
end

subaxis(2,2,4)
plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k'); hold on;
title(extractBefore(param_idx(9,1),'_'))
axis off
for i = 1:size(Tg2_runs,2)
    temp_run = Tg2_runs{1,i};
    temp_circ = Tg2_circ(i,1);
    ismember(unique_circ,Tg2_circ(1,1));
    plot(temp_run(:,1),temp_run(:,2),'Color',cm_magma(ismember(unique_circ,Tg2_circ(i,1)),:),'LineWidth',3); hold on;
end

subaxis(2,2,1)
plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k'); hold on;
title(extractBefore(param_idx(45,1),'_'))
axis off
for i = 1:size(Wt1_runs,2)
    temp_run = Wt1_runs{1,i};
    temp_circ = Wt1_circ(i,1);
    ismember(unique_circ,Wt1_circ(1,1));
    plot(temp_run(:,1),temp_run(:,2),'Color',cm_magma(ismember(unique_circ,Wt1_circ(i,1)),:),'LineWidth',3); hold on;
end

subaxis(2,2,3)
plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k'); hold on;
title(extractBefore(param_idx(29,1),'_'))
axis off
for i = 1:size(Wt2_runs,2)
    temp_run = Wt2_runs{1,i};
    temp_circ = Wt2_circ(i,1);
    ismember(unique_circ,Wt2_circ(1,1));
    plot(temp_run(:,1),temp_run(:,2),'Color',cm_magma(ismember(unique_circ,Wt2_circ(i,1)),:),'LineWidth',3); hold on;
end

ax = axes;
c = colorbar(ax);
ax.Visible = 'off';
colormap(cm_magma)

%%
param_idx=params.subID;
fig=figure;
fig.Color = [1 1 1];
for i=1:2:48
    [f,x] = ecdf(params.stop_cue_proximity{i,1}');
    if contains(param_idx(i),'Tg')
        plot(x,f,'LineWidth',3,'Color','#601a4a');
    else
        plot(x,f,'LineWidth',3,'Color','#9c9eb5');
    end
    hold on;
end
grid on
set(gca,'FontSize',14,'FontWeight','bold')

ylabel('Proportion')
xlabel('Circuity')


%% Stops by stop duration (size) and time in trial (color)
for i=1:size(params.tsStop,1)
    for ii=1:size(params.tsStop{i},2)
        params.tsStopIdx{i}(1,ii)=params.tsStop{i}{1,ii}(1,1);
    end
end


tg1=vertcat(params.stopCenter{contains(param_idx,'Tg') & contains(param_idx,'D1')});
tg1_c=horzcat(params.timeStopped{contains(param_idx,'Tg') & contains(param_idx,'D1')});
tg1_ts=horzcat(params.tsStopIdx{contains(param_idx,'Tg') & contains(param_idx,'D1')});

wt1=vertcat(params.stopCenter{contains(param_idx,'WT') & contains(param_idx,'D1')});
wt1_c=horzcat(params.timeStopped{contains(param_idx,'WT') & contains(param_idx,'D1')});
wt1_ts=horzcat(params.tsStopIdx{contains(param_idx,'WT') & contains(param_idx,'D1')});

fig=figure;
fig.Color=[1 1 1];
subaxis(1,2,2)
plot(sin(0:2*pi/1000:2*pi)*102,cos(0:2*pi/1000:2*pi)*102,'k'); hold on;
scatter(tg1(:,1),tg1(:,2),cell2mat(tg1_c)',cell2mat(tg1_c)','Filled')
title('Tg')
axis off
colormap(magma(255))
caxis([min([cell2mat(tg1_c),cell2mat(wt1_c)]),max([cell2mat(tg1_c),cell2mat(wt1_c)])])

subaxis(1,2,1)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
scatter(wt1(:,1),wt1(:,2),cell2mat(wt1_c)',cell2mat(wt1_c)','Filled')
title('Wt')
axis off
colorbar
colormap(magma(255))
caxis([min([cell2mat(tg1_c),cell2mat(wt1_c)]),max([cell2mat(tg1_c),cell2mat(wt1_c)])])

param_idx=params.subID;
fig=figure;
fig.Color = [1 1 1];
for i=1:2:48
    [f,x] = ecdf(params.timeStopped{i,1}');
    if contains(param_idx(i),'Tg')
        plot(x,f,'LineWidth',3,'Color','#601a4a');
    else
        plot(x,f,'LineWidth',3,'Color','#9c9eb5');
    end
    hold on;
end
grid on
set(gca,'FontSize',14,'FontWeight','bold')

ylabel('Proportion')
xlabel('Stop Duration (s)')

%% Home Base Figures

for i=1:size(params.HBcenter,1)
    for ii = 1:size(params.hbOcc{i},2)
        move_slow{i,1}(1,ii) = cell2mat(params.HBclass{i}(1,ii)) > .75;
    end
end

%num of HB that meet above
for i = 1 : length (params.hbOcc)
    numHB(i,1)= sum(move_slow{i});
end

for i=1:size(params.HBcenter,1)
    [~,occIdx(i,1)]=max(params.hbOcc{i}(move_slow{i})); %max given moving slow most of time;
end

wt_d1 = contains(param_idx,{'WT'}) & contains(param_idx,{'D1'});
tg_d1 = contains(param_idx,{'Tg'}) & contains(param_idx,{'D1'});

wt_hb_idx = occIdx(wt_d1,1);
tg_hb_idx = occIdx(tg_d1,1);

% convert hexidecimal color to rgb for fill
str = '#9c9eb5';
wt_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#601a4a';
tg_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

fig = figure;
fig.Color = [1 1 1];
subaxis(1,2,1)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
title('Wt')
test_Wt = params(wt_d1,:);
for i = 1:size(wt_hb_idx,1)
    %     temp = test_Wt.HBBound{i}{1,wt_hb_idx(i,1)};
    for ii = 1:size(test_Wt.HBBound{i},2)
        temp = test_Wt.HBBound{i}{1,ii};
        fill(temp(:,1),temp(:,2),wt_color,'FaceAlpha',.5,'EdgeAlpha',.5,'EdgeColor',wt_color)
    end
end
axis image
axis off

subaxis(1,2,2)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
title('Tg')
test_Tg = params(tg_d1,:);
for i = 1:size(wt_hb_idx,1)
    %     temp = test_Tg.HBBound{i}{1,tg_hb_idx(i,1)};
    for ii = 1:size(test_Tg.HBBound{i},2)
        temp = test_Tg.HBBound{i}{1,ii};
        fill(temp(:,1),temp(:,2),tg_color,'FaceAlpha',.5,'EdgeAlpha',.5,'EdgeColor',tg_color)
    end
end
axis image
axis off

%% Figure 5 probe test

% At least one visit to cue location (proportion of animals)
fig = figure;
fig.Color = [1 1 1];
x = categorical({'Wt','Tg'});
x = reordercats(x,{'Wt','Tg'});
b = bar(x,[prop_cue_day2_Wt/12,prop_cue_day2_Tg/12],'EdgeAlpha',0);
b.FaceColor = 'flat'
b.CData(1,:) = wt_color;
b.CData(2,:) = tg_color;
box off
ylabel('Proportion')
xlabel('Group')
set(gca,'FontSize',14,'FontWeight','bold')

% Example Trajectories to Cue location

wt_d2 = contains(param_idx,{'WT'}) & contains(param_idx,{'D2'});
tg_d2 = contains(param_idx,{'Tg'}) & contains(param_idx,{'D2'});

fig = figure;
fig.Color = [1 1 1];
subaxis(1,2,1)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
title('Wt')

test_Wt = params(wt_d2,:);
for i = 1:size(test_Wt,1)
    if isnan(test_Wt.tsEntry_cue{i})
        continue
    end
    temp = test_Wt.in_cue_zonw_idx{i};
    first = find(temp, 1, 'first');
    x = test_Wt.backCM{i}(1:first,1);
    y = test_Wt.backCM{i}(1:first,2);
    plot(x,y,'LineWidth',2)
    plot(x(1,1),y(1,1),'*r'); hold on
end
axis image
axis off

subaxis(1,2,2)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
title('Tg')

test_Tg = params(tg_d2,:);
for i = 1:size(test_Tg,1)
    if isnan(test_Tg.tsEntry_cue{i})
        continue
    end
    temp = test_Tg.in_cue_zonw_idx{i};
    first = find(temp, 1, 'first');
    x = test_Tg.backCM{i}(1:first,1);
    y = test_Tg.backCM{i}(1:first,2);
    plot(x,y,'LineWidth',2)
    plot(x(1,1),y(1,1),'*r'); hold on
end
axis image
axis off

%% Old plot versions

%% Surface plots for occupancy (sum to get total time and divide by 60 to get into minutes

% create occupancy map for each animal
for i = 1:length(params.subID)
    x = params.backCM{i}(:,1);
    y = params.backCM{i}(:,2);
    [params.occMap{i},map] = OF.occ_map(x,y,202,3,30);
end

tg1=sum(cat(3,params.occMap{contains(param_idx,'Tg') & contains(param_idx,'D1')}),3)/12;
wt1=sum(cat(3,params.occMap{contains(param_idx,'WT') & contains(param_idx,'D1')}),3)/12;
tg2=sum(cat(3,params.occMap{contains(param_idx,'Tg') & contains(param_idx,'D2')}),3)/12;
wt2=sum(cat(3,params.occMap{contains(param_idx,'WT') & contains(param_idx,'D2')}),3)/12;

%Use meshgrid to serve as basis for logical mask.
imageSizeX = size(tg1,1);
imageSizeY = size(tg1,2);
[columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

clear imageSizeX imageSizeY

% Next create the circle in the image.
centerX = median(1:size(tg1,1)); centerY = median(1:size(tg1,2)); radius = median(1:size(tg1,2));
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

clear rowsInImage columnsInImage

tg1(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
wt1(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
tg2(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
wt2(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN

imAlpha1=ones(size(tg1));
imAlpha1(isnan(tg1))=0;
imAlpha2=ones(size(wt1));
imAlpha2(isnan(wt1))=0;
imAlpha3=ones(size(tg2));
imAlpha3(isnan(tg2))=0;
imAlpha4=ones(size(wt2));
imAlpha4(isnan(wt2))=0;

fig=figure; fig.Color=[1 1 1];
subplot(1,2,2)
ax = imagesc(tg1,'AlphaData',imAlpha1); colormap(magma(255));axis xy;hold on; box off; axis image;shading flat;
caxis([min(min(min(tg1,wt1))),max(max(max(tg1,wt1)))])
set(gca, 'CLim', [0, 20])
title('Tg')
colorbar('southoutside');
axis off

subplot(1,2,1)
imagesc(wt1,'AlphaData',imAlpha2); colormap(magma(255));axis xy;hold on; box off; axis image;shading flat;
caxis([min(min(min(tg1,wt1))),max(max(max(tg1,wt1)))])
set(gca, 'CLim', [0, 20])
title('Wt')
axis off

