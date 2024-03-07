close all
clearvars;

addpath('./BedformIdentification');
addpath('./BedformSeparation');

% Load an example transect
load example_transect.mat

%%%%% Separate bathymetric data  representing multiscale bedforms %%%%%%%%%
% Set parameters: 
d_x = 15;               % The half span of the LOESS smoother. Approximately 3-5 times the secondary bedform length
window = 10;            % Window in which exact break location can be set. Around 0.2 times the primary dune length
cutoff_slope = 0.05;    % Cutoff slope, default 0.03

% Separate the bed elevation profile. z_separation is the resulting fitted
% line, idx_sigmoid indicates where the fitted line is determined by a
% sigmoid fit. This information is used in determining the maximum lee side
% angle of primary dunes. s is the longitudinal coordinate along the
% curvilinear grid.
[z_separation, idx_sigmoid]= bedformSeparation(s,z,d_x,window,cutoff_slope);

% For larger datasets, the bedform separation can be run in parallel, to
% decrease computational time. 

% Plot
figure; 
subplot(2,1,1); hold on; plot(s,z,'k'); 
    plot(s,z_separation,'r');
    xlim([min(s) max(s)]); ylabel('z [m]'); xlabel('s [m]');
subplot(2,1,2); hold on; plot(s, z-z_separation','k'); % Secondary bedforms
    plot(s,zeros(length(z),1),'r')
    xlim([min(s) max(s)]); ylabel('z [m]'); xlabel('s [m]');


%%%%% Identify bedforms and calculate properties %%%%%%%%%%%%%%%%%%%%%%%%%%
% Put data in the correct structure: 
data.grids = s; % Longitudinal coordinate along curvilinear grid. dimensions: [number of datapoints along the profile, number of transects]
data.gridn = n; % Cross-sectional coordinate. 
data.gridx = x; % In case of curvilinear grid, used to compute the length along grid
data.gridy = y; % In case of curvilinear grid, used to compute the length along grid 
data.gridz = z;% dimensions: [number of datapoints along the profile, number of transects, number of timesteps]
data.gridz_separation = z_separation'; % dimensions: [number of datapoints along the profile, number of transects, number of timesteps]
data.idx_sigmoid = idx_sigmoid'; % Required for determining max primary lee slope

timesteps = 1;
transects = 1;
propsSec = propsSecondary(data,timesteps, transects);
propsPrim = propsPrimary(data, timesteps, transects);

% Filtering parameters can be adjusted in the script.
remove = 'yes'; % yes or no: determines whether filtered bedforms are removed or only flagged. 
[propsSec, propsPrim, perc] = propsFilter(propsSec, propsPrim, data, timesteps, transects, remove);

% Plot profile with identified troughs and tops. 
transect=1;time=1;
lw=1;
color1 = [27,158,119]./255; color2 = [217,95,2]./255; color3 = [117,112,179]./255;
marker_sec    = 'o'; markersize_sec = 2;
marker_prim   = 'o'; markersize_prim = 5;

figure; hold on;
plot(data.grids(transect,:),data.gridz(transect,:,time),'LineWidth',lw, 'color','k');
plot(data.grids(transect,:),data.gridz_separation(transect,:,time),'LineWidth',lw,'color',color1);
plot(data.grids(transect, propsSec(transect,time).cr),...
    data.gridz(transect, propsSec(transect,time).cr,time),'.',...
    'marker',marker_sec,'markersize',markersize_sec, 'MarkerFaceColor',color2, 'color', color2);
plot(data.grids(transect, propsSec(transect,time).tr1),...
    data.gridz(transect, propsSec(transect,time).tr1,time),'.',...
    'marker',marker_sec,'markersize',markersize_sec, 'MarkerFaceColor',[1,1,1], 'color', color2);
plot(data.grids(transect, propsSec(transect,time).tr2),...
    data.gridz(transect, propsSec(transect,time).tr2,time),'.',...
    'marker',marker_sec,'markersize',markersize_sec, 'MarkerFaceColor',[1,1,1], 'color', color2);

plot(data.grids(transect,propsPrim(transect,time).cr),...
    data.gridz_separation(transect,propsPrim(transect,time).cr,time),'.',...
    'marker',marker_prim,'markersize',markersize_prim, 'MarkerFaceColor',color3,'color', color3);
plot(data.grids(transect,propsPrim(transect,time).tr1),...
    data.gridz_separation(transect,propsPrim(transect,time).tr1,time),'.',...
    'marker',marker_prim,'markersize',markersize_prim, 'MarkerFaceColor',[1,1,1],'color', color3);
plot(data.grids(transect,propsPrim(transect,time).tr2),...
    data.gridz_separation(transect,propsPrim(transect,time).tr2,time),'.',...
    'marker',marker_prim,'markersize',markersize_prim, 'MarkerFaceColor',[1,1,1],'color', color3);

