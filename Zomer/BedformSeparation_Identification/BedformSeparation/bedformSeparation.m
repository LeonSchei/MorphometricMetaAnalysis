function [zSeparation,idxSigmoid] = bedformSeparation(x,z,dx,window,cutoff_slope)
% x is a vector with regularly spaced x-coordinates and should be positive in
% upstream direction. 
%
% z is a vector containing the z-coordinates, corresponding to x. z is positive in upward
% direction. 
%
% dx is the half-span of the smoother. A higher values of dx results in a smoother LOESS
% curve. dx is approximately 3-5 times the secondary dune length
%
% cutoff_slope is the slope of the primary lee side slope at which a break
% needs to be introduced. Default 0.03 m/m
%
% window is approximately 0.1-0.2 times the primary dune length
%
%  Judith Zomer, 2021-11-28

% input x z should be vector x * 1. 
if size(x,2)>1;x=x';end
if size(z,2)>1;z=z';end

% Remove nan-values from edges. 
zSeparation = nan(length(z),1);
idxSigmoid=nan(length(z),1);
idxnotnan = ~isnan(z);
z=z(idxnotnan); x=x(idxnotnan);

%%%% STEP 1: Unmodified LOESS curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unmodified_loess_curve = loess_interp(x,z,0,x,dx,2);%

%%%% STEP 2: Approximate break locations based on unmodified LOESS curve %%
[troughs_z,troughs_x] = findpeaks(-unmodified_loess_curve,x);troughs_z=-troughs_z;
[~,crests_x] = findpeaks(unmodified_loess_curve,x);

break_trough=[];break_crest=[];
z_slope = gradient(unmodified_loess_curve)./gradient(x); % Local slope of z


for i=1:size(troughs_x)
 lss = find(crests_x > troughs_x(i)); % Find for each trough the upstream crest
 if ~isempty(lss) 
 lss = lss(1); 
    if max(z_slope((x >= troughs_x(i)) & (x <= crests_x(lss))))>cutoff_slope  % If max lee slope is higher than cutoff_slope, breakpoint at trough location
     break_trough(i)=troughs_x(i);
     break_crest(i) = crests_x(lss); 
    end
 end
end

% If there are no breakpoints, apply LOESS, else, continue
if isempty(break_trough) 
    z_separation = loess_interp(x,z,0,x,dx,2);
    idx_sigmoid  = zeros(length(x),1);
else; break_crest(break_crest==0)=[];break_trough(break_trough==0)=[]; 
    
%%%% STEP 3: Update break location based on unfiltered bed elevation profile %%%%%%%%
% Find minimum z (unfiltered data) within specified distance ('window') upstream from 
% breakpoint
for i=1:length(break_trough)
  idx = (x>break_trough(i) & x<break_trough(i)+window);
  x_subset = x(idx); z_subset=z(idx);
  x_subset = x_subset(z_subset == min(z_subset));
  break_trough(i) = x_subset(end); % breakpoints is updated
end
         
%%%% Move breakpoint upstream if slope z is very low %%%%%%%%%%%%%%
for i=1:length(break_trough)
    idx = (x >= break_trough(i) & x < break_trough(i)+window);
    x_subset = x(idx); z_subset=z(idx);
    for ii=2:length(x_subset)
        if (z_subset(ii)-z(x==break_trough(i)))/(x_subset(ii)-break_trough(i)) < cutoff_slope 
            break_trough(i) = x_subset(ii);
        end
    end
end
clear x_subset z_subset
    
%%%% Move breakpoint crest downstream if slope z is very low %%%%%%%%%%%%%%
for i=1:length(break_crest)
    idx = (x > break_crest(i)-window & x < break_crest(i));
    x_subset = x(idx); z_subset=z(idx);
    for ii=2:length(x_subset)
        if abs((z_subset(ii)-z(x==break_crest(i)))/(x_subset(ii)-break_crest(i))) < cutoff_slope  
            break_crest(i) = x_subset(ii);
        end
    end
end
clear x_subset z_subset
    
%%%% STEP 4, 5 :Break up the transect based on x-values in breakpoints and 
%%%% apply LOESS stepwise. At steep lee side slopes, a sigmoid function is
%%%% fitted
    idx=find(break_crest-break_trough <=0); % Remove break if trough is further upstream than crest.
    break_crest(idx)=[];break_trough(idx)=[];

    break_trough = [min(x), break_trough, max(x)]; % Add 0 and max x to breakpoints
    [break_trough,idx] = unique(break_trough); 
    
    break_crest = [min(x), break_crest,max(x)]; % Add 0 and max x to breakpoints
    break_crest = break_crest(idx);
    
    x_fitted=[];z_fitted=[];z_separation=[];idx_sigmoid=[];
    for i=1:length(break_trough)-1
      % Break to crest    
      % Fit sigmoid function to dune lee side. 
      if ~(i==1)                                                            
      x_tofit = x(x >= break_trough(i) & x <= break_crest(i));             % Select x,z data to be fitted
      z_tofit = z(x >= break_trough(i) & x <= break_crest(i));
      sigmoidfun = @(b,x)(b(1) + b(2)./(1 + exp(-b(3).*(x-b(4)))));         % define sigmoid function
      Fsumsquares = @(b) ...
          sum((sigmoidfun(b,x_tofit(1:end)) - z_tofit(1:end)).^2);          % Set the objective funtion as sum of squares of the residuals
      options = optimoptions(@fmincon, 'Algorithm','interior-point',...
          'MaxFunctionEvaluations',3000,'Display','off');
      beta0=[z_tofit(1) z_tofit(end)-z_tofit(1) 0.5 x_tofit(1)];            % Initial values beta
      [beta] = fmincon(Fsumsquares, beta0, [],[],[],[],...                  % Find optimum
          [beta0(1)-4*beta0(2) beta0(2)*0.5 0 beta0(4)-8*window],...        % lower bounds
          [beta0(1)+beta0(2)/4 beta0(2)*4   3 beta0(4)+4*window],...
          [],options);   
      z_separation_subset=sigmoidfun(beta,x_tofit);                         % Calculate the fitted lines based on found beta and the sigmoid function
     
      % Connect sigmoid function to previous loess (stoss) and next (crest)
      z_separation_slope=gradient(z_separation_subset)./gradient(x_tofit);
      maxslope_idx = find(z_separation_slope==max(z_separation_slope));
      lowslope_lee2crest_idx = find(z_separation_slope < ...
          max([0.5*max(z_separation_slope),cutoff_slope]));
      lowslope_lee2crest_idx(lowslope_lee2crest_idx<maxslope_idx(1))=[];                % only keep the part near the crest (idx) 
      lowslope_lee2crest_idx=min(lowslope_lee2crest_idx);                               % idlowslope=upper part of sigmoid that is removed. 
     
      x_tofit_lee2crest = x_tofit(lowslope_lee2crest_idx:end);                          % upper part of sigmoid that is removed.
      z_tofit_lee2crest = z_tofit(lowslope_lee2crest_idx:end);                          % upper part of sigmoid that is removed.
      
      x_tofit(lowslope_lee2crest_idx:end)=[];
      z_separation_subset(lowslope_lee2crest_idx:end)=[];
      z_tofit(lowslope_lee2crest_idx:end)=[];                                           % Upper part of sigmoid removed from x,z, and lpass
      
      if ~isempty(x_tofit)
        lowslope_trough2lee_idx = find(z_separation_slope<...
            max([0.5*max(z_separation_slope),cutoff_slope]));
        lowslope_trough2lee_idx(lowslope_trough2lee_idx>maxslope_idx(1))=[];
        lowslope_trough2lee_idx=max(lowslope_trough2lee_idx);                           %find lower part of sigmoid that should be removed. 

        if isempty(lowslope_trough2lee_idx)
            lowslope_trough2lee_idx=1;
        end
        lowslope_trough2lee_idx = min(lowslope_trough2lee_idx, length(x_tofit)); 
        
        % Update idlowslope if z value is lower then absolute trough
        % location.
        zmintrough = min(z(x >= break_trough(i)-window & x <= break_trough(i)+window));
        id_zmintrough = find(z_separation_subset < zmintrough);
        if ~isempty(id_zmintrough); lowslope_trough2lee_idx=max(lowslope_trough2lee_idx, id_zmintrough(end));end

        % x,z of trough where LOESS has to be fitted to connect to the sigmoid fit 
        xx_idtrough = x_tofit(1:lowslope_trough2lee_idx);
        zz_idtrough = z_tofit(1:lowslope_trough2lee_idx); 
        
        % Remove values from x_fitted and z_separation if lower then z mintrough
        % starting from last index. 
        if lowslope_trough2lee_idx < 5
                    xx_idtrough=[x_fitted(end-(5-lowslope_trough2lee_idx):end); xx_idtrough];
                    zz_idtrough=[z_fitted(end-(5-lowslope_trough2lee_idx):end); zz_idtrough];
                    z_separation(end-(5-lowslope_trough2lee_idx):end)=[];
                    x_fitted(end-(5-lowslope_trough2lee_idx):end)=[];
                    z_fitted(end-(5-lowslope_trough2lee_idx):end)=[];
                    idx_sigmoid(end-(5-lowslope_trough2lee_idx):end)=[];
        end
        while 1 
            if length(z_separation)>1        
                if z_separation(end)<zmintrough
                    xx_idtrough=[x_fitted(end);xx_idtrough];
                    zz_idtrough=[z_fitted(end);zz_idtrough];
                    z_separation(end)=[];x_fitted(end)=[];z_fitted(end)=[];idx_sigmoid(end)=[]; 
                else
                    break; 
                end
            else
                break;
            end
        end
        
        % LOESS fit section between downstream LOESS curve and sigmoid,
        % force trough.. 
       xxtrough = [repmat(x_fitted(end),100,1); xx_idtrough;repmat(x_tofit(lowslope_trough2lee_idx),100,1)];% ;
       zztrough = [repmat(z_separation(end),100,1); zz_idtrough;repmat(z_separation_subset(lowslope_trough2lee_idx),100,1)];
        
       lpasstrough = loess_interp(xxtrough,zztrough,0,xxtrough,dx,2);% LOess fit between previous and sigmoid function; here nans?

       [~,idx] = unique(xxtrough);idx=[find(xxtrough==0);idx];
       xxtrough = xxtrough(idx);
       zztrough = zztrough(idx);
       lpasstrough = lpasstrough(idx);
      
       x_fitted = [x_fitted;xxtrough;x_tofit(lowslope_trough2lee_idx:end)]; % connect fit to loess_trough and sigmoid function. 
       z_fitted = [z_fitted;zztrough;z_tofit(lowslope_trough2lee_idx:end)];
       z_separation = [z_separation;lpasstrough;z_separation_subset(lowslope_trough2lee_idx:end)];
       idx_sigmoid = [idx_sigmoid; zeros(length(xxtrough),1);ones(length(x_tofit(lowslope_trough2lee_idx:end)),1)];
      end
      end
      % Crest to next break
      if i==1 % Fit first section of the bed elevation profile, up to the first trough
        x_tofit = [x(x >= break_trough(i) & x <= break_trough(i+1))];
        z_tofit = [z(x >= break_trough(i) & x <= break_trough(i+1))]; 
      elseif isempty(x_tofit) % If the lee side slope was not steep enough, apply LOESS from trough to next trough
          x_tofit = [x(x >= break_trough(i) & x <= break_trough(i+1))];
          z_tofit = [z(x >= break_trough(i) & x <= break_trough(i+1))]; 
          x_tofit = [repmat(x_fitted(end),100,1);x_tofit];% ;
          z_tofit = [repmat(z_separation(end),100,1);z_tofit];% );
      else % If sigmoid function was fitted, apply LOESS from crest to next trough
         x_tofit = [x(x >= break_crest(i) & x <= break_trough(i+1))];
         z_tofit = [z(x >= break_crest(i) & x <= break_trough(i+1))];
         x_tofit = [repmat(x_fitted(end),100,1); x_tofit_lee2crest;x_tofit];% ;
         z_tofit = [repmat(z_separation(end),100,1); z_tofit_lee2crest;z_tofit];% );
      end
      z_separation_subset = loess_interp(x_tofit,z_tofit,0,x_tofit,dx,2);
      [~,idx] = unique(x_tofit);idx=[find(x_tofit==0);idx];
      x_tofit = x_tofit(idx);
      z_tofit = z_tofit(idx);
      z_separation_subset = z_separation_subset(idx);
      x_fitted = [x_fitted;x_tofit];
      z_fitted = [z_fitted;z_tofit];
      z_separation = [z_separation;z_separation_subset];
      idx_sigmoid = [idx_sigmoid; zeros(length(x_tofit),1)];
    end
    [~, ~, subs] = unique(x_fitted);
    z_separation = accumarray(subs, z_separation, [], @mean);
    idx_sigmoid = accumarray(subs, idx_sigmoid, [], @min);
end
zSeparation(idxnotnan)=z_separation;
idxSigmoid(idxnotnan)=idx_sigmoid;
end


