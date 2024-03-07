function [propsSec, propsPrim, perc] = propsFilter(propsSec, propsPrim, data, timesteps, transects, remove)
for time = timesteps
for i=transects
    propsPrim(i,time).filter = [];
    propsPrim(i,time).filter = propsPrim(i,time).H < 0.25 | propsPrim(i,time).H>4 | propsPrim(i,time).L > 200 | ...
    propsPrim(i,time).L<25 | propsPrim(i,time).A > 0.2 | propsPrim(i,time).A < 0.005 | ...
    propsPrim(i,time).maxLee<0.05;
   
    perc.primary(i,time) = sum(propsPrim(i,time).filter)./length(propsPrim(i,time).H).*100;
    
    if isequal(remove, 'yes')
        propsPrim(i,time).tr1   = propsPrim(i,time).tr1(~propsPrim(i,time).filter);    
        propsPrim(i,time).tr2   = propsPrim(i,time).tr2(~propsPrim(i,time).filter);
        propsPrim(i,time).cr    = propsPrim(i,time).cr(~propsPrim(i,time).filter);

        propsPrim(i,time).H = propsPrim(i,time).H(~propsPrim(i,time).filter);
        propsPrim(i,time).fH = propsPrim(i,time).fH(~propsPrim(i,time).filter);
        propsPrim(i,time).L = propsPrim(i,time).L(~propsPrim(i,time).filter);
        propsPrim(i,time).Lee = propsPrim(i,time).Lee(~propsPrim(i,time).filter,:);
        propsPrim(i,time).fLee = propsPrim(i,time).fLee(~propsPrim(i,time).filter,:);
        propsPrim(i,time).Stoss = propsPrim(i,time).Stoss(~propsPrim(i,time).filter,:);
        propsPrim(i,time).fStoss = propsPrim(i,time).fStoss(~propsPrim(i,time).filter,:);
        propsPrim(i,time).Depth = propsPrim(i,time).Depth(~propsPrim(i,time).filter); 
        propsPrim(i,time).A = propsPrim(i,time).A(~propsPrim(i,time).filter);
        propsPrim(i,time).fA = propsPrim(i,time).fA(~propsPrim(i,time).filter); 
       % propsPrim(i,time).Time = propsPrim(i,time).Time(~propsPrim(i,time).filter);
        clear propsPrim(i,time).filter;
    end


    propsSec(i,time).filter = propsSec(i,time).H < 0.05 | propsSec(i,time).H>0.75 | propsSec(i,time).L > 25 | ...
        propsSec(i,time).L<0.5 | propsSec(i,time).A > 0.2 | propsSec(i,time).A < 0.005 | ...
        data.gridz(i,propsSec(i,time).cr,time)'<data.gridz(i,propsSec(i,time).tr1,time)'-0.01 | ...
        data.gridz(i,propsSec(i,time).cr,time)'<data.gridz(i,propsSec(i,time).tr2,time)'-0.01 | ...
        propsSec(i,time).maxLee<0.05;

     perc.secondary(i,time) = sum(propsSec(i,time).filter)./length(propsSec(i,time).H).*100; 
 
    if isequal(remove, 'yes')
        propsSec(i,time).tr1 = propsSec(i,time).tr1(~propsSec(i,time).filter);% Contains IDs
        propsSec(i,time).tr2 = propsSec(i,time).tr2(~propsSec(i,time).filter);
        propsSec(i,time).cr = propsSec(i,time).cr(~propsSec(i,time).filter);

        propsSec(i,time).H = propsSec(i,time).H(~propsSec(i,time).filter);
        propsSec(i,time).fH = propsSec(i,time).fH(~propsSec(i,time).filter);
        propsSec(i,time).L = propsSec(i,time).L(~propsSec(i,time).filter);
        propsSec(i,time).Lee = propsSec(i,time).Lee(~propsSec(i,time).filter,:);
        propsSec(i,time).fLee = propsSec(i,time).fLee(~propsSec(i,time).filter,:);
        propsSec(i,time).Stoss = propsSec(i,time).Stoss(~propsSec(i,time).filter,:);
        propsSec(i,time).fStoss = propsSec(i,time).fStoss(~propsSec(i,time).filter,:);
        propsSec(i,time).Depth = propsSec(i,time).Depth(~propsSec(i,time).filter);
        propsSec(i,time).A = propsSec(i,time).A(~propsSec(i,time).filter);
        propsSec(i,time).fA = propsSec(i,time).fA(~propsSec(i,time).filter); 
      %  propsSec(i,time).Time = propsSec(i,time).Time(~propsSec(i,time).filter);
        clear propsSec(i,time).filter;
    end 
    
end
end
end