function [troughs, crests] = findExtremes(y, z, upcross, downcross)
    crossings = [upcross; downcross];
    if isempty(crossings); troughs=[]; crests=[]; else % If crossings is empty
    crossings = sortrows(crossings,1);

    % Find crests and troughs
    crests = []; troughs = [];
    for i = 1:length(crossings(:,1))-1
        if crossings(i+1,3)-crossings(i,3) == 1 % Then top
            f = []; f = find((y > crossings(i,1)) & (y < crossings(i+1,1)));
            dist2f = y(f); z_f = z(f); provo0f = z(f);
            idx_crest = find(z_f==max(z_f)); idx_crest=idx_crest(1);
            z_crest = z_f(idx_crest);
            ycrest = dist2f(idx_crest);
            crests = [crests; ycrest z_crest 3]; % 

        elseif crossings(i+1,3)-crossings(i,3) == -1 % Then trough
            f = find((y > crossings(i,1)) & (y < crossings(i+1,1)));
            dist2f = y(f); z_f = z(f); provo0f = z(f);
            z_troughtemp = min(z_f);
            idx_trough = find(z_f==min(z_f)); idx_trough=idx_trough(1);
            z_trough = z_f(idx_trough);
            ytrough = dist2f(idx_trough);      
            troughs = [troughs; ytrough z_trough 4]; %
        end
    end

    % Check if crests and troughs alternate
    if ~isempty(crests) && ~isempty(troughs)
        if abs(diff([size(crests,1),size(troughs,1)]))<=1
        else; error('The amount of troughs does not agree with the amount of crests!')
        end
    end
    extremes = [crests; troughs];

    if ~isempty(extremes)
        extremes = sortrows(extremes,1);
        diffs = abs(diff(extremes(:,3)));
        if ~isempty(diffs~=1) 
        else; error('Troughs and crests do not alternate')      
        end
    end
    end
end


