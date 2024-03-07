function [upcross, downcross] = findCrossings(y, z)

% Find zero-crossings
    upcross = []; downcross = [];
    for i = 1:(length(z)-1)
        if z(i) * z(i+1) < 0 % zero-crossing if smaller than zero
            fit = polyfit(y(i:i+1),z(i:i+1),1);
            crossing = -fit(2)/fit(1); % x-coordinaat van de up- of downcrossing
            if fit(1) > 0 % rico positief, dan upcrossing, krijgt label 1 mee
                upcross = [upcross; crossing 0 1];
            elseif fit(1) < 0 % rico negatief, dan downcrossing, krijgt label 2 mee
                downcross = [downcross; crossing 0 2];
            end
        end
    end

% Check if up- and downcrossings alternate

    crossings = [upcross; downcross];
    if isempty(crossings); upcross=[]; downcross=[]; else % Allow for no crossings
    crossings = sortrows(crossings,1);

    if ~isempty(upcross) && ~isempty(downcross)
        if abs(diff([size(upcross,1),size(downcross,1)]))<=1
        else; error('The amount of upcrossings does not agree with the amount of downcrossings!')
        end
    end
    
    diffs = abs(diff(crossings(:,3)));
    if ~isempty(diffs~=1) 
    else; error('Up- and downcrossings do not alternate')      
    end 

    upcross = crossings(crossings(:,3)==1,:);
    downcross = crossings(crossings(:,3)==2,:);
    end 
end


