% this is an example of how to process bathymetric data with the method
% presented in Lefebvre et al. (2021)
% it can easily be adapted to rapidly process many files, for example with
% a for-loop

% example data, with many thanks to Wasserstraßen- und Schifffahrtsamt 
% Weser-Jade-Nordsee, Standort Bremerhaven for allowing us to share the data
load('Weser_20090210_km47-55')   % not yet available
  
% xr and yr should be (measures of) longitude and latitude
% if you have x as streamwise and y as crosswise, use the following line
%xr=yr; yr=xr; zr=zr;

% It is also assumed that size(xr,1)<size(xr,2), 
if size(xr,1)<size(xr,2)
    xr=xr';    yr=yr';    zr=zr';
end

% getting the crestlines and troughlines
MC=get_crest_troughlines(xr,yr,zr,'');  % threshold left as blank

% %%%%%%%% uncomment below to create a figure with crest and throughlines 
figure
pcolor(xr,yr,zr)
axis equal
shading flat
% % drawing crestlines
for n=1:MC.CL.NumObjects
    line(xr(MC.CL.PixelIdxList{n}),yr(MC.CL.PixelIdxList{n}),'color','k','linewidth',2)
    %line(xr(MC.CL.PixelIdxList{n}),MC.CL.smoothed{n},'color','r','linewidth',2)
end
% drawing troughlines
for n=1:MC.TL.NumObjects
    line(xr(MC.TL.PixelIdxList{n}),yr(MC.TL.PixelIdxList{n}),'color','w','linewidth',2)
end
% drawing bifurcations
bifurc=find(MC.pos_bif==1);
line(xr(bifurc),yr(bifurc),'markerfacecolor','r','marker','o','linestyle','none')

% getting the slip faces
steep_face=tidal_steep_faces(xr,yr,zr,15,45);

% getting properties
BP=tidal_bedform_properties(xr,yr,zr,MC.crtr,steep_face);
