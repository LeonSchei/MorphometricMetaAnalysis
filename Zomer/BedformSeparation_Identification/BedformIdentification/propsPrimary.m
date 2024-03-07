function propsPrim = propsPrimary(data, timesteps, transects)

propsPrim(1:size(data.gridz,1),1:size(data.gridz,3))=struct('H',[],'fH',[],...
    'L',[],'Depth',[],'tr1',[],'cr',[],'tr2',[],'Lee',[],'fLee',[],'maxLee',[],'maxfLee',[],'fStoss',[],'Stoss',[],'A',[],'fA',[]);

for time = timesteps
    for i = transects
        x = data.grids(i,:)';
        z = data.gridz_separation(i,:,time)';
        z = interp1(x(~isnan(z)),z(~isnan(z)),x,'linear','extrap');
        %gridz_separation = movmean(data.gridz_tmean(i,:),1000)';gridz_separation=gridz_separation+(mean(z)-mean(gridz_separation));
        %gridz_separation=loess_interp(x,z,0,x,data.primarylength*12,2);
        gridz_separation = movmean(z,4000);gridz_separation=gridz_separation+(mean(z)-mean(gridz_separation));
        
        [upcross, downcross] = findCrossings(x, z-gridz_separation);
        [troughs, crests]    = findExtremes(x, z, upcross, downcross);
        
% Focus on trough to crest height, so all crest locations downstream of
% first trough are removed % And remove crest location upstream of trough (dune should be trough-crest-trough)    
    crests(crests(:,1) < min(troughs(:,1)),:)=[];
    crests(crests(:,1) > max(troughs(:,1)),:)=[];
        
    crest.z = crests(:,2); crest.s = crests(:,1);[~,crest.id] = ismember(crest.s, x);
    trough.z = troughs(:,2); trough.s = troughs(:,1);[~,trough.id] = ismember(trough.s, x);
    
 % Height is crest_z - average height_z
    fHeight = crest.z-nanmean([trough.z(1:end-1),trough.z(2:end)],2);
    Height = data.gridz(i,crest.id,time)'-nanmean([data.gridz_separation(i,trough.id(1:end-1),time)',data.gridz_separation(i,trough.id(2:end),time)'],2);
 
 Wrong = Height < 0.25 | fHeight < 0.25;
 Wrong_x = crest.s(Wrong);
  
 % Delete upcrossings and downcrossings next to Wrong_x
 % First mark for all, then delete, because might overlap(?)--> no don't
 % overlap. 
 upcross_up = upcross;
 downcross_up = downcross;
 for id=1:length(Wrong_x)
     % Upcross: lower then Wrong_x(i), then closest. 
 upcross_dif = Wrong_x(id) - upcross_up(:,1);
 upcross_dif(upcross_dif<0)=nan;idx = find(upcross_dif==min(upcross_dif));
 upcross_up(idx,:)=[];
 
 downcross_dif = Wrong_x(id) - downcross_up(:,1);
 downcross_dif(downcross_dif>0)=nan;idx = find(downcross_dif==-min(abs(downcross_dif)));
 downcross_up(idx,:)=[];
 end
 clear Height fHeight heightx Wrong Wrong_x SNct SNtc
 
[troughs, crests]                = findExtremes(x, z, upcross_up, downcross_up);% should make this function myself-does not find the 'absolute' extremes...       
          
% Focus on trough to crest height, so all crest locations downstream of
% first trough are removed % And remove crest location upstream of trough (dune should be trough-crest-trough)    
    crests(crests(:,1) < min(troughs(:,1)),:)=[];
    crests(crests(:,1) > max(troughs(:,1)),:)=[];
        
    crest.z = crests(:,2); crest.s = crests(:,1);[~,crest.id] = ismember(crest.s, x);
    trough.z = troughs(:,2); trough.s = troughs(:,1);[~,trough.id] = ismember(trough.s, x);
    
 % Height is crest_z - average height_z
    fHeight = crest.z-nanmean([trough.z(1:end-1),trough.z(2:end)],2);
    Height = data.gridz(i,crest.id,time)'-nanmean([data.gridz_separation(i,trough.id(1:end-1),time)',data.gridz_separation(i,trough.id(2:end),time)'],2);
  
    
% SLOPE (Should I average around grid nodes?) CHECK
    fSlopet = diff(z)./...
        sqrt((data.gridx(i,2:end)'-data.gridx(i,1:end-1)').^2+...
        (data.gridy(i,2:end)'-data.gridy(i,1:end-1)').^2);
    Slopet  = diff(data.gridz_separation(i,:,time)')./...
        sqrt((data.gridx(i,2:end)'-data.gridx(i,1:end-1)').^2+...
        (data.gridy(i,2:end)'-data.gridy(i,1:end-1)').^2);
    
    fSlope(1)=fSlopet(1); fSlope(2:length(fSlopet))=nanmean([fSlopet(1:end-1),fSlopet(2:end)],2);
    fSlope(length(fSlopet)+1)=fSlopet(end);
    Slope(1)=Slopet(1); Slope(2:length(Slopet))=nanmean([Slopet(1:end-1),Slopet(2:end)],2);
    Slope(length(Slopet)+1)=Slopet(end);   
    
     
    % Length: use dx.
    dx = sqrt((data.gridx(i,2:end)-data.gridx(i,1:end-1)).^2+...
        (data.gridy(i,2:end)-data.gridy(i,1:end-1)).^2);
    Length = nan(size(fHeight));
    for idx=1:length(trough.id)-1
        Length(idx)=sum(dx(trough.id(idx):trough.id(idx+1)));
    end
    
    % Depth: should be weighted average... 
    Depth = nan(size(fHeight));
    for idx=1:length(trough.id)-1
       Depth(idx) = sum(data.gridz(i,trough.id(idx):trough.id(idx+1)).*dx(trough.id(idx):trough.id(idx+1)))./Length(idx);
    end
    
    Lee = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        Lee(idx,1) = nanmean(Slope(trough.id(idx):crest.id(idx)));
    end
     
    fLee = nan(size(fHeight,1),1);
    for idx=1:length(trough.id)-1
        fLee(idx,1) = nanmean(fSlope(trough.id(idx):crest.id(idx)));
    end    

    maxLee=nan(size(fHeight,1),1);    
    maxfLee=nan(size(fHeight,1),1);
    for idx = 1:length(trough.id)-1
        if sum(data.idx_sigmoid(i,trough.id(idx):crest.id(idx)),'omitnan')==0
            maxLee(idx) = max(Slope(trough.id(idx):crest.id(idx)),[],'omitnan');
            maxfLee(idx) = max(fSlope(trough.id(idx):crest.id(idx)),[],'omitnan');
        else
            Slope_subset = Slope(trough.id(idx):crest.id(idx));idx_sigmoid_subset = data.idx_sigmoid(i,trough.id(idx):crest.id(idx));
            fSlope_subset = fSlope(trough.id(idx):crest.id(idx));
            maxLee(idx) = max(Slope_subset(idx_sigmoid_subset==1),[],'omitnan');
            maxfLee(idx) = max(fSlope_subset(idx_sigmoid_subset==1),[],'omitnan');
        end
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
    
propsPrim(i,time).tr1       = trough.id(1:end-1);% Contains IDs
propsPrim(i,time).tr2       = trough.id(2:end);
propsPrim(i,time).cr        = crest.id;

propsPrim(i,time).H         = Height;
propsPrim(i,time).fH        = fHeight;
propsPrim(i,time).L         = Length;
propsPrim(i,time).Lee       = Lee;
propsPrim(i,time).fLee      = fLee;
propsPrim(i,time).maxLee    = maxLee;
propsPrim(i,time).maxfLee   = maxfLee;
propsPrim(i,time).Stoss     = Stoss;
propsPrim(i,time).fStoss    = fStoss;
propsPrim(i,time).Depth     = Depth;
propsPrim(i,time).A         = A;
propsPrim(i,time).fA        = fA; 
%propsPrim(i,time).Time(:,1) = data.time(i,crest.id,time);
    end
end
end