function [Dune] = IdentifyDunes4(lx,lz)
%------------------------------------------------------------------
%	Find location and dimensions for Dunes belonging to 
%   length classes t according to Ashley et al. (1990)
%   along an arbitrary transect defined by x and z vectors
%
%   t = 1 : Ripples          (0 m < L < 0.6 m)
%   t = 2 : Small Dunes      (0.6 m < L < 5 m)
%   t = 3 : Medium Dunes     (5 m > L > 10 m)
%   t = 4 : Large Dunes      (10 m < L 100 m)
%   t = 5 : Very Large Dunes (100 m < L < oo)
%
%   This version includes the loop over Ashley's classes and was
%   complemented by a double-detection based on the function 'unique'.
%
%   Detailed information about this approach is given in:
%   Scheiber, Leon, et al. - "Robust methods for the decomposition 
%   and interpretation of compound dunes applied to a complex hydro-
%   morphological setting." Earth Surface Processes and Landforms,
%   46.2 (2021): 478-489.
%
%------------------------------------------------------------------
%% Set global definitions
DL = [0.01 0.6   5 10 100 500];    % D(t)Length (m) / small / medium / large / very large
DH = 0.0677.*DL.^0.8098;           % D(t)Height (m) / small / medium / large / very large

p = 0.15;   % Relative height boundary (default)
Vacc = 0.05; % Vertical MBES accuracy = lower detection limit
dx = lx(2)-lx(1);   % Horizontal resolution (m)

%% Pre-Allocate the following D(t) variables
h = [];       % Height (Relative)
a = [];       % Asymmetry
l = [];       % Length

pp = [];      % Peak-Position
ci = [];      % Peak-Index
pv = [];      % Peak-Value (Absolute)

tv = [];      % Trough-Value
tp = [];      % Trough-Index
ti = [];      % Trough-Position

for t = 1:5
    
    %% Find Minima and Maxima
    pvn = []; ppn = [];
    [pv  pp ] = findpeaks( lz,lx','MinPeakProminence',p*DH(t),'Annotate','extents');
    [pvn ppn] = findpeaks(-lz,lx','MinPeakProminence',p*DH(t),'Annotate','extents');
    
    %% Find Peakpositions in original Z-vector
    id = []; idn = [];
    for r = 1:numel(pp)                 % Loop over all found peaks
        id(r) =  find(lx == pp(r));     % Find x and y value on original Z
        pv(r) =  lz(id(r));
    end
    for rn = 1:numel(ppn)                 % Loop over all found peaks
        idn(rn) =  find(lx == ppn(rn));   % Find x and y value on original Z
        pvn(rn) =  lz(idn(rn));
    end
    
    %% Find vallies
    lv = []; rv = []; d = [];
    for N = 1:numel(id) % loop through peaks
        
        rxl = idn(idn<id(N)); % ANGEPASST ppn statt idn
        rxr = idn(idn>id(N));
        dl = min(id(N) - rxl); % ANGEPASST
        dr = min(rxr - id(N));
        dl(dl*dx>=DL(t+1)) = []; % ANGEPASST: relevant für dx < 1 m (...)
        dr(dr*dx>=DL(t+1)) = [];
        
        if isempty(dl) || isempty(dr)
            lv(N) = 0;
            rv(N) = 0;
        else
            [~, lv(N)] = min(abs(lx(id(N)-dl) - lx.'));
            [~, rv(N)] = min(abs(lx(dr+id(N)) - lx.'));
        end
    end
    
    %% Delete fails    
    out = find(lv==0 | rv==0);
    lv(out) = []; rv(out) = [];
    pp(out) = []; pv(out) = []; 
    id(out) = []; ci = id;
    
    %% Create trough vectors
    ti = sort([lv, rv]); % combine left and right valley into vector
    tp = lx(ti);         % define absolute trough positions
    tv = lz(ti);         % extract valley heights (m) from originial Z
    
    %% Asymmetries according to Zorndt et al., 2011 | eq. (1)
    a = [];
    for k = 1:numel(pp)
        grad_d =  []; % gradient variable
        grad_d = gradient(lz(ti(k):ti(k+1)),lx(ti(k):ti(k+1))); %
        a(k) = (numel(grad_d(grad_d<=0)) - numel(grad_d(grad_d>0))) / numel(grad_d); % Asymetry
    end
    
    %% Calculate D(t) heights and lengths
    h1 = []; h2 = []; h3 = []; l1 = []; l2 = []; l3 = []; alpha = []; % Ex = []; Ey = [];
    for k = 1:numel(pv)
        kk = 2*k-1;
        % Mean value
        h3(k) = abs(pv(k) - (tv(kk)+tv(kk+1))/2);
        
        % OR Vertical distance: Peak - interpolated trend line
        h(k) = abs(pv(k) - (tv(kk)+ (tv(kk+1)-tv(kk)) * (pp(k)-tp(kk)) / (tp(kk+1)-tp(kk))) );
        
        % OR Orthogonal distance according to Wesseling & Wilbers (2000)
        Ax = tp(kk);    Ay = tv(kk);
        Bx = pp(k);     By = pv(k);
        Cx = tp(kk+1);  Cy = tv(kk+1);
        Dx = Bx;
        if Cy > Ay
            Dy = Ay + (Cy-Ay) * ( (Bx-Ax) / (Cx-Ax) );
            h1(k) = By - Dy;
            alpha(k) = atan( (Cx-Ax) / (Cy-Ay) );
            h2(k) = h1(k) * (sin(alpha(k)));
        else
            Dy = Cy + (Ay-Cy) * ( (Cx-Bx) / (Cx-Ax) );
            h1(k) = By - Dy;
            alpha(k) = atan( (Cx-Ax) / (Ay-Cy) );
            h2(k) = h1(k) * (sin(alpha(k)));
        end
             
    % Transect baseline:
    l2(k) = sqrt( (tp(kk+1) - tp(kk))^2 + ...
                 (tv(kk+1) - tv(kk))^2 );     
    % Sum of horiz. distances:
    l1(k) = sum(diff(lx(ti(kk):ti(kk+1))));   
    % Sum of baseline increments:
    l3(k) = sqrt(sum(diff(lx(ti(kk):ti(kk+1))))^2 + ... 
                    (tv(kk+1) - tv(kk))^2 );
    end
    
    %% Allow for instrument accuracy
    out = find(h1 <= Vacc | l1 <= 5*dx);
    h1(out) = []; h2(out) = []; h3(out) = [];
    l1(out) = []; l2(out) = []; l3(out) = []; 
    a(out) = []; alpha(out) = [];
    id(out) = []; pv(out) = []; pp(out) = []; 
    ti([out*2-1 out*2]) = []; tv([out*2-1 out*2]) = []; tp([out*2-1 out*2]) = [];  
        
    %% Create struct
    D(t).H1 = h1;                 % Height (Relative)
    D(t).H2 = h2;                 % Height (Relative)
    D(t).H3 = h3;                 % Height (Relative)
    D(t).L1 = l1;                 % Length
    D(t).L2 = l2;                 % Length
    D(t).L3 = l3;                 % Length
    D(t).A = a;                   % Asymmetry
    D(t).D = pv - 0.5.*h1;       % Mean Depth
    D(t).alpha = alpha;           % Baseline inclination
    
    D(t).ci = id;                 % Peak-Index
    D(t).ti1 = ti(1:2:end-1);     % First Trough-Index
    D(t).ti2 = ti(2:2:end);       % Second Trough-Index
    
    D(t).lx = lx;                 % Distance-Vector
    D(t).lz = lz;                 % Height-Vector
    
end

%% Extract all unique value pairs
ind = find(arrayfun(@(D) ~isempty(D.H1),D));
Du = unique([[D(ind).H1]', [D(ind).H2]', [D(ind).H3]', [D(ind).L1]', [D(ind).L2]', [D(ind).L3]',...
            [D(ind).A]', [D(ind).D]', [D(ind).ci]', [D(ind).ti1]', [D(ind).ti2]', [D(ind).alpha]'],'rows');

%% Combine all type results
if ~isempty(Du)
    Dune.H = Du(:,3)';     % Average dune height (peak to trendline)
    Dune.H1 = Du(:,1)';     % Vertical dune height (peak to trendline)
    Dune.H2 = Du(:,2)';     % Orthogonal dune height (peak to trendline)
    Dune.H3 = Du(:,3)';     % Average dune height (peak to up- and downstream trough)
    Dune.L = Du(:,4)';     % Dune length (1)
    Dune.L1 = Du(:,4)';     % Dune length (1)
    Dune.L2 = Du(:,5)';     % Dune length (1)
    Dune.L3 = Du(:,6)';     % Dune length (1)
    Dune.A = Du(:,7)';      % Dune asymmetry
    Dune.D = Du(:,8)';      % Mean dune depth
    Dune.alpha = Du(:,12)'; % Baseline inclination
    
    Dune.ti1 = Du(:,10)';   % Left trough indices
    Dune.ci = Du(:,9)';     % Crest indices
    Dune.ti2 = Du(:,11)';   % Right trough indices
    
    Dune.lx = D(1).lx;      % Profile distance
    Dune.lz = D(1).lz;      % Profile height
else
    return
end

end
