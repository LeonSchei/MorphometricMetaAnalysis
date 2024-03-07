function propsSec = propsSecondary(data,timesteps, transects)
% Inititialize
propsSec(1:size(data.gridz,1),1:size(data.gridz,3))=struct('H',[],'fH',[],...
    'L',[],'Depth',[],'tr1',[],'cr',[],'tr2',[],'Lee',[],'fLee',[],'maxLee',[],'maxfLee',[],'fStoss',[],...
    'Stoss',[],'A',[],'fA',[],'Time',[]);

for time = timesteps
    for i = transects
        x = data.grids(i,:)';
        z = data.gridz(i,:,time)' - data.gridz_separation(i,:,time)';
        zu = data.gridz(i,:,time)';        
        z = interp1(x(~isnan(z)),z(~isnan(z)),x,'linear','extrap');%remove nan values
        
% Identify troughs and crests based on zero-crossing.
    [upcross, downcross] = findCrossings(x, z);%
    [troughs, crests]    = findExtremes(x, z, upcross, downcross);
  
% All crest locations downstream of the first trough are removed 
% And remove crest location upstream of trough (dune should be trough-crest-trough)    
    crests(crests(:,1) < min(troughs(:,1)),:)=[];
    crests(crests(:,1) > max(troughs(:,1)),:)=[];
        
% z, z_unfiltered, s and id in crest, trough structure
    crest.z = crests(:,2); crest.s = crests(:,1);[~,crest.id] = ismember(crest.s, x);crest.zu=zu(crest.id);
    trough.z = troughs(:,2); trough.s = troughs(:,1);[~,trough.id] = ismember(trough.s, x);trough.zu=zu(trough.id);
  
 % Plot transect with crest and trough locations
 %   figure; plot(x,zu,'k');hold on; plot(troughs(:,1),zu(trough.id),'.r')
 %   plot(crest.s,zu(crest.id),'.b'); plot(x,data.gridz_separation(i,:,time)','color',[0.5 0.5 0.5])
    
 % Height is crest_z - average height_z
   fHeight = crest.z-nanmean([trough.z(1:end-1),trough.z(2:end)],2);
   Height = data.gridz(i,crest.id,time)'-nanmean([data.gridz(i,trough.id(1:end-1),time)',data.gridz(i,trough.id(2:end),time)'],2);

  
% Remove crossings that lead to a bedform height smaller than 0.02 m.
   Wrong=Height<0.02 | fHeight < 0.02;%
   Wrong_x = crest.s(Wrong);
  
 % Delete upcrossings and downcrossings next to Wrong_x
  upcross_idfalse = nan(length(upcross),1);
  downcross_idfalse = nan(length(downcross),1);
 
 for id=1:length(Wrong_x)
     % Upcross: lower then Wrong_x(i), then closest. 
 upcross_dif = Wrong_x(id) - upcross(:,1);
 upcross_dif(upcross_dif<0)=nan;idx = find(upcross_dif==min(upcross_dif));
 upcross_idfalse(idx)=1;
 
 downcross_dif = Wrong_x(id) - downcross(:,1);
 downcross_dif(downcross_dif>0)=nan;idx = find(downcross_dif==-min(abs(downcross_dif)));
 downcross_idfalse(idx)=1;
 end
 
 upcross_up = upcross(upcross_idfalse~=1,:);
 downcross_up = downcross(downcross_idfalse~=1,:);
 
clear Height fHeight heightx Wrong Wrong_x SNct SNtc
 
% Find troughs and crests after filtering crossings
[troughs, crests]                = findExtremes(x, z, upcross_up, downcross_up);
        
% First troughs are removed 
% And remove crest location upstream of trough (dune should be trough-crest-trough)    
    crests(crests(:,1) < min(troughs(:,1)),:)=[];
    crests(crests(:,1) > max(troughs(:,1)),:)=[];
        
% z, z_unfiltered, s and id in crest, trough structure
    crest.z = crests(:,2); crest.s = crests(:,1);[~,crest.id] = ismember(crest.s, x);crest.zu=zu(crest.id);
    trough.z = troughs(:,2); trough.s = troughs(:,1);[~,trough.id] = ismember(trough.s, x);trough.zu=zu(trough.id);
     
    % Height is crest_z - average height_z
    fHeight = crest.z-nanmean([trough.z(1:end-1),trough.z(2:end)],2);
    Height = data.gridz(i,crest.id,time)'-nanmean([data.gridz(i,trough.id(1:end-1),time)',data.gridz(i,trough.id(2:end),time)'],2);
 
 
% Slope
    fSlopet = diff(z)./...
        sqrt((data.gridx(i,2:end)'-data.gridx(i,1:end-1)').^2+...
        (data.gridy(i,2:end)'-data.gridy(i,1:end-1)').^2);
    Slopet  = diff(data.gridz(i,:,time)')./...
        sqrt((data.gridx(i,2:end)'-data.gridx(i,1:end-1)').^2+...
        (data.gridy(i,2:end)'-data.gridy(i,1:end-1)').^2);
    
    fSlope(1)=fSlopet(1); fSlope(2:length(fSlopet))=nanmean([fSlopet(1:end-1),fSlopet(2:end)],2);
    fSlope(length(fSlopet)+1)=fSlopet(end);
    Slope(1)=Slopet(1); Slope(2:length(Slopet))=nanmean([Slopet(1:end-1),Slopet(2:end)],2);
    Slope(length(Slopet)+1)=Slopet(end);   
        
% Length: use dx in case of curvilinear grid.
    dx = sqrt((data.gridx(i,2:end)-data.gridx(i,1:end-1)).^2+...
        (data.gridy(i,2:end)-data.gridy(i,1:end-1)).^2);
    Length = nan(size(fHeight));
    for idx=1:length(trough.id)-1
        Length(idx)=sum(dx(trough.id(idx):trough.id(idx+1)));
    end
    
    % Depth
    Depth = nan(size(fHeight));
    for idx=1:length(trough.id)-1
       Depth(idx) = sum(data.gridz(i,trough.id(idx):trough.id(idx+1)).*dx(trough.id(idx):trough.id(idx+1)))./Length(idx);
    end
    
    % Lee slope
    Lee = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        Lee(idx,1) = nanmean(Slope(trough.id(idx):crest.id(idx)));
    end
     
    fLee = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        fLee(idx,1) = nanmean(fSlope(trough.id(idx):crest.id(idx)));
    end
    
    maxLee = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        maxLee(idx,1) = max(Slope(trough.id(idx):crest.id(idx)),[],'omitnan');
    end
     
    maxfLee = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        maxfLee(idx,1) = max(fSlope(trough.id(idx):crest.id(idx)),[],'omitnan');
    end 
    
    
    % Stoss
    Stoss = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        Stoss(idx,1) = nanmean(Slope(crest.id(idx):trough.id(idx+1)));
    end
     
    fStoss = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        fStoss(idx,1) = nanmean(fSlope(crest.id(idx):trough.id(idx+1)));
    end   
    
A = Height./Length;
fA = fHeight./Length; 
    
propsSec(i,time).tr1 = trough.id(1:end-1);% Contains IDs
propsSec(i,time).tr2 = trough.id(2:end);
propsSec(i,time).cr  = crest.id;

propsSec(i,time).H        = Height;
propsSec(i,time).fH       = fHeight;
propsSec(i,time).L        = Length;
propsSec(i,time).Lee      = Lee;
propsSec(i,time).fLee     = fLee;
propsSec(i,time).maxLee   = maxLee;
propsSec(i,time).maxfLee  = maxfLee;
propsSec(i,time).Stoss    = Stoss;
propsSec(i,time).fStoss   = fStoss;
propsSec(i,time).Depth    = Depth;
propsSec(i,time).A        = A;
propsSec(i,time).fA       = fA; 
%propsSec(i,time).Time(:,1)= data.time(i,crest.id,time);
    end
end