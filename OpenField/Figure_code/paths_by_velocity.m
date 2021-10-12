%% Example paths for 12 months small open field across days
save_path = '/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/figs/MatlabFigs/small_of_paths_by_vel/';

% load df with small path coords and velocity
load('/Users/lauraberkowitz/github/TgF344-AD_Open_Field/notebooks/data/df_small_OF.mat')
% subset the df to retain only 12 month paths 
df_sub = df([df.time_point{:}] == 12,:); 

%% create colormap the size of the unique values in data -
% now each line can be colorcoded on the same scale by the velocity value
temp_vel = [];
for path = 1:length(df_sub.rat)
    temp_vel = [temp_vel; df_sub.velocity{path}];
end

temp_vel(temp_vel < 0) = 0;

unique_vel = unique([temp_vel]','sorted');
cm_magma = colormap(magma(size(unique_vel,2)));

% loop through rats and plot 3 days

for rat = 1:3:length(df_sub.rat)
    
    fig = figure;
    fig.Color = [1 1 1];
    
    [fig_width_in, fig_height_in] = set_size('thesis', .75, [1,1]);
    set(fig,'Position',[835 270 fig_width_in fig_height_in])
    
    subplot(1,3,1)
    % maze boundary
    plot(sin(0:pi/360:2*pi)*(76.5/2),cos(0:pi/360:2*pi)*(76.5/2),'k'); hold on;
    title('Day 1')
    axis off
    
    % plot path colorcoded by velocity
    plot(df_sub.x_cm{rat},df_sub.y_cm{rat},'Color',[.75 .75 .75],'LineWidth',2); hold on;
    scatter(df_sub.x_cm{rat},df_sub.y_cm{rat},50,df_sub.velocity{rat},'filled');
    axis image
    ax = axes;
    c = colorbar(ax);
    ax.Visible = 'off';
    
    subplot(1,3,2)
    % maze boundary
    plot(sin(0:pi/360:2*pi)*(76.5/2),cos(0:pi/360:2*pi)*(76.5/2),'k'); hold on;
    title('Day 2')
    axis off
    
    % plot path colorcoded by velocity
    plot(df_sub.x_cm{rat+1},df_sub.y_cm{rat+1},'Color',[.75 .75 .75],'LineWidth',2); hold on;
    scatter(df_sub.x_cm{rat+1},df_sub.y_cm{rat+1},50,df_sub.velocity{rat+1},'filled');
    axis image
    ax = axes;
    c = colorbar(ax);
    ax.Visible = 'off';
    
    subplot(1,3,3)
    % maze boundary
    plot(sin(0:pi/360:2*pi)*(76.5/2),cos(0:pi/360:2*pi)*(76.5/2),'k'); hold on;
    title('Day 3')
    axis off
    
    % plot path colorcoded by velocity
    plot(df_sub.x_cm{rat+2},df_sub.y_cm{rat+2},'Color',[.75 .75 .75],'LineWidth',2); hold on;
    scatter(df_sub.x_cm{rat+2},df_sub.y_cm{rat+2},50,df_sub.velocity{rat+2},'filled');
    axis image
    ax = axes;
%     c = colorbar(ax);
    ax.Visible = 'off';
    colormap(cm_magma)
    
    saveas(fig,[save_path,filesep,num2str(df_sub.rat(rat)),'_12mo_paths_velocity','.png'],'png')
    
    close all
    
end



