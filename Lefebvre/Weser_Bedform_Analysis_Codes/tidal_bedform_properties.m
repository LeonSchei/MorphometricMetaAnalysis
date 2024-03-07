function BP=tidal_bedform_properties(xr,yr,zr,crtr,steep_face)
%BP=tidal_bedform_properties(xr,yr,zr,crtr,steep_face)
% function to analyse tidal bedform properties
% input are coordinates along the waterway, matrix with the crest and
% trough position (possibly coming MC.crtr from get_crest_troughlines.m) and
% the position of the steep face (coming from  tidal_steep_faces.m)
%
% output is a structure with all the properties
% 
% Alice Lefebvre, 2021 alefebvre@marum.de

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, create matrix with the same size as input where the properties
% will be stored
Li=nan(size(xr));                       % bedform length from trough to trough
Hi=nan(size(xr));                       % bedform height from mean(trough) to crest
Lebblee=nan(size(xr));                  % length ebb lee side
Hebblee=nan(size(xr));                  % height ebb lee side
MSebblee=nan(size(xr));                 % mean angle ebb lee side
MSmaxebblee=nan(size(xr));              % max angle ebb lee side
xposMSmaxebblee=nan(size(xr));          % horizontal distance between crest and position max angle ebb lee side
zposMSmaxebblee=nan(size(xr));          % vertical distance between crest and position max angle ebb lee side (after Cisneros et al 2020)
Lfloodlee=nan(size(xr));                % length flood lee side
Hfloodlee=nan(size(xr));                % height flood lee side
MSfloodlee=nan(size(xr));               % mean angle flood lee side
MSmaxfloodlee=nan(size(xr));            % max angle flood lee side
xposMSmaxfloodlee=nan(size(xr));        % horizontal distance between crest and position max angle flood lee side
zposMSmaxfloodlee=nan(size(xr));        % vertical distance between crest and position max angle ebb flood side (after Cisneros et al 2020)
xcrest=nan(size(xr));                   % x pos crest
ycrest=nan(size(xr));                   % y pos cerst
zcrest=nan(size(xr));                   % depth crest
xtrough1=nan(size(xr));                 % x pos trough 1
ytrough1=nan(size(xr));                 % y pos trough 1
ztrough1=nan(size(xr));                 % z pos trough 1
xtrough2=nan(size(xr));                 % x pos trough 2
ytrough2=nan(size(xr));                 % y pos trough 2
ztrough2=nan(size(xr));                 % z pos trough 2
meanbedformdepth=nan(size(xr));         % mean bedform depth
MESP=nan(size(xr));                     % mean slope ebb steep face
HESP=nan(size(xr));                     % mean height ebb steep face
LESP=nan(size(xr));                     % mean length ebb steep face
x_esp_beg=nan(size(xr));                % beginning ebb steep face
x_esp_end=nan(size(xr));                % end ebb steep face
x_esp_begN=nan(size(xr));               % normalised pos beginning ebb steep face
x_esp_endN=nan(size(xr));               % normalised pos end ebb steep face
MFSP=nan(size(xr));                     % mean slope flood steep face
HFSP=nan(size(xr));                     % mean height flood steep face
LFSP=nan(size(xr));                     % mean length flood steep face
x_fsp_beg=nan(size(xr));                % beginning flood steep face
x_fsp_end=nan(size(xr));                % end flood steep face
x_fsp_begN=nan(size(xr));               % normalised pos beginning flood steep face
x_fsp_endN=nan(size(xr));               % normalised pos end flood steep face
isthereSP=nan(size(xr));                % is there a steep face 0 = no steep face, 1 = ebb steep face, -1 = flood steep face and 2 = both steep face

% looking per transect where crests and troughs are. If it is found:
% trough - crest - trough, then  properties are calcualted and stored at 
% the position of the crest
for N=1:size(xr,2)
    CT=find(crtr(:,N)~=0);
    for i=1:length(CT)-2
        
        tr1_pos=CT(i);
        cr_pos=CT(i+1);
        tr2_pos=CT(i+2);
        
        % bedform only if find trough crest trough 
        if crtr(tr1_pos,N)==-1  && crtr(cr_pos,N)==1 && crtr(tr2_pos,N)==-1
            
            fn=find(isnan(zr(tr1_pos:tr2_pos,N))==1);
            
            if length(fn)<3 && tr2_pos>tr1_pos+3        % only if less than 3 NaNs in the bathymetry and at least 4 points
                zr(tr1_pos:tr2_pos,N) = fillmissing(zr(tr1_pos:tr2_pos,N),'linear');
                steep_face.slopemag(tr1_pos:tr2_pos,N) = fillmissing(steep_face.slopemag(tr1_pos:tr2_pos,N),'linear');
                steep_face.zr_sw(tr1_pos:tr2_pos,N) = fillmissing(steep_face.zr_sw(tr1_pos:tr2_pos,N),'linear');
                %cr_pos=find(zr(tr1_pos:tr2_pos,N)==max(zr(cr_pos-5:cr_pos+5,N)),1)+tr1_pos-1;
  
                Li(CT(i+1),N)=sqrt((xr(tr2_pos,N)-xr(tr1_pos,N))^2+(yr(tr2_pos,N)-yr(tr1_pos,N))^2);
                Hi(CT(i+1),N)=zr(cr_pos,N)-(mean([zr(tr1_pos,N) zr(tr2_pos,N)]));
                Lebblee(CT(i+1),N)=sqrt((xr(tr2_pos,N)-xr(cr_pos,N))^2+(yr(tr2_pos,N)-yr(cr_pos,N))^2);
                Hebblee(CT(i+1),N)=zr(cr_pos,N)-zr(tr2_pos,N);
                Lfloodlee(CT(i+1),N)=sqrt((xr(cr_pos,N)-xr(tr1_pos,N))^2+(yr(cr_pos,N)-yr(tr1_pos,N))^2);
                Hfloodlee(CT(i+1),N)=-zr(tr1_pos,N)+zr(cr_pos,N);
                MSfloodlee(CT(i+1),N)=mean(steep_face.slopemag(tr1_pos:cr_pos-1,N));
                MSebblee(CT(i+1),N)=mean(steep_face.slopemag(cr_pos:tr2_pos-1,N));
                
                MSmaxfloodlee(CT(i+1),N)=max(steep_face.slopemag(tr1_pos:cr_pos-1,N));
                pmf=find(steep_face.slopemag(tr1_pos:cr_pos-1,N)==MSmaxfloodlee(CT(i+1),N));
                zposMSmaxfloodlee(CT(i+1),N)=(steep_face.zr_sw(tr1_pos+pmf-1,N)-zr(tr1_pos,N))/Hfloodlee(CT(i+1),N);
                xposMSmaxfloodlee(CT(i+1),N)=sqrt((steep_face.xr_sw(tr1_pos+pmf-1,N)-xr(cr_pos,N))^2+(steep_face.yr_sw(tr1_pos+pmf-1,N)-yr(cr_pos,N))^2)/Lfloodlee(CT(i+1),N);
                
                MSmaxebblee(CT(i+1),N)=min(steep_face.slopemag(cr_pos:tr2_pos-1,N));
                pme=find(steep_face.slopemag(cr_pos:tr2_pos-1,N)==MSmaxebblee(CT(i+1),N))-1;
                zposMSmaxebblee(CT(i+1),N)=(steep_face.zr_sw(cr_pos+pme,N)-zr(tr2_pos,N))/Hebblee(CT(i+1),N);
                xposMSmaxebblee(CT(i+1),N)=sqrt((steep_face.xr_sw(cr_pos+pme,N)-xr(cr_pos,N))^2+(steep_face.yr_sw(cr_pos+pme,N)-yr(cr_pos,N))^2)/Lebblee(CT(i+1),N);
                
                xcrest(CT(i+1),N)=xr(cr_pos,N);
                ycrest(CT(i+1),N)=yr(cr_pos,N);
                zcrest(CT(i+1),N)=zr(cr_pos,N);
                xtrough1(CT(i+1),N)=xr(tr1_pos,N);
                ytrough1(CT(i+1),N)=yr(tr1_pos,N);
                ztrough1(CT(i+1),N)=zr(tr1_pos,N);
                xtrough2(CT(i+1),N)=xr(tr2_pos,N);
                ytrough2(CT(i+1),N)=yr(tr2_pos,N);
                ztrough2(CT(i+1),N)=zr(tr2_pos,N);
                meanbedformdepth(CT(i+1),N)=nanmean(zr(tr1_pos:tr2_pos,N));
                %                 figure
                %                 plot(xr(tr1_pos-5:tr2_pos+5,N),zr(tr1_pos-5:tr2_pos+5,N),'.k-')
                %                 hold on
                %                 plot(xr(tr1_pos,N),zr(tr1_pos,N),'og')
                %                 plot(xr(tr2_pos,N),zr(tr2_pos,N),'og')
                %                 plot(xr(cr_pos,N),zr(cr_pos,N),'*r')
                %                 plot(steep_face.xr_sw(tr1_pos-5:tr2_pos+5,N),steep_face.zr_sw(tr1_pos-5:tr2_pos+5,N),'xk')
                %                 plot(steep_face.xr_sw(tr1_pos+pmf-1,N),steep_face.zr_sw(tr1_pos+pmf-1,N),'*m')
                %                 plot(steep_face.xr_sw(cr_pos+pme,N),steep_face.zr_sw(cr_pos+pme,N),'*g')
                
                % finding ebb steep face
                ESP=find(steep_face.ebb(tr1_pos:tr2_pos,N)>0);
                % calculating ebb steep face properties
                if isempty(ESP)==0
                    
                    n1=ESP(1)+tr1_pos-1;
                    n2=ESP(end)+tr1_pos;
                    
                    MESP(CT(i+1),N)=mean(steep_face.slopemag(ESP+tr1_pos-1,N));
                    HESP(CT(i+1),N)=zr(n1,N)-zr(n2,N);
                    x_esp_beg(CT(i+1),N)=sqrt((xr(n1,N)-xr(cr_pos,N))^2+(yr(n1,N)-yr(cr_pos,N))^2);
                    x_esp_end(CT(i+1),N)=sqrt((xr(n2,N)-xr(cr_pos,N))^2+(yr(n2,N)-yr(cr_pos,N))^2);
                    LESP(CT(i+1),N)=x_esp_end(CT(i+1),N)-x_esp_beg(CT(i+1),N);
                    x_esp_begN(CT(i+1),N)=x_esp_beg(CT(i+1),N)./Lebblee(CT(i+1),N);
                    x_esp_endN(CT(i+1),N)=x_esp_end(CT(i+1),N)./Lebblee(CT(i+1),N);
                    
                    isthereSP(CT(i+1),N)=1;
                    %                     plot(steep_face.xr_sw(ESP(1)+tr1_pos-1:ESP(end)+tr1_pos-1,N),steep_face.zr_sw(ESP(1)+tr1_pos-1:ESP(end)+tr1_pos-1,N),'x-r')
                    %                     plot(xr(n1,N),zr(n1,N),'dr')
                    %                     plot(xr(n2,N),zr(n2,N),'dr')
                clear n1 n2
                end
                
                % finding flood steep face
                FSP=find(steep_face.flood(tr1_pos:tr2_pos,N)>0);
                % calculating ebb steep face properties
                if isempty(FSP)==0
                    n1=FSP(end)+tr1_pos;
                    n2=FSP(1)+tr1_pos-1;
                    
                    MFSP(CT(i+1),N)=mean(steep_face.slopemag(tr1_pos+FSP-1,N));
                    HFSP(CT(i+1),N)=zr(n1,N)-zr(n2,N);
                    x_fsp_beg(CT(i+1),N)=sqrt((xr(n1,N)-xr(cr_pos,N))^2+(yr(n1,N)-yr(cr_pos,N))^2);
                    x_fsp_end(CT(i+1),N)=sqrt((xr(n2,N)-xr(cr_pos,N))^2+(yr(n2,N)-yr(cr_pos,N))^2);
                    LFSP(CT(i+1),N)=x_fsp_end(CT(i+1),N)-x_fsp_beg(CT(i+1),N);
                    x_fsp_begN(CT(i+1),N)=x_fsp_beg(CT(i+1),N)./Lfloodlee(CT(i+1),N);
                    x_fsp_endN(CT(i+1),N)=x_fsp_end(CT(i+1),N)./Lfloodlee(CT(i+1),N);
                    
                    isthereSP(CT(i+1),N)=-1;
                    %                     plot(steep_face.xr_sw(FSP(1)+tr1_pos-1:FSP(end)+tr1_pos-1,N),steep_face.zr_sw(FSP(1)+tr1_pos-1:FSP(end)+tr1_pos-1,N),'x-r')
                    %                     plot(xr(n1,N),zr(n1,N),'dr')
                    %                     plot(xr(n2,N),zr(n2,N),'dr')
                end
                
                if isempty(FSP)==0 && isempty(ESP)==0
                    isthereSP(CT(i+1),N)=2;
                end
                
                if isempty(FSP) && isempty(ESP)
                    isthereSP(CT(i+1),N)=0;
                end
                
            end
        end
    end
end
  
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BP.L=Li;                                % bedform length from trough to trough
BP.H=Hi;                                % bedform height from mean(trough) to crest
BP.Lebblee=Lebblee;                     % length ebb lee side        
BP.Hebblee=Hebblee;                     % height ebb lee side
BP.MSebblee=MSebblee;                   % mean angle ebb lee side
BP.MSmaxebblee=MSmaxebblee;             % max angle ebb lee side
BP.xposMSmaxebblee=xposMSmaxebblee;     % position max angle ebb lee side (after Cisneros et al 2020)
BP.zposMSmaxebblee=zposMSmaxebblee;     % position max angle ebb lee side (after Cisneros et al 2020)
BP.Lfloodlee=Lfloodlee;                 % length flood lee side
BP.Hfloodlee=Hfloodlee;                 % height flood lee side
BP.MSfloodlee=MSfloodlee;               % mean angle flood lee side
BP.MSmaxfloodlee=MSmaxfloodlee;         % max angle flood lee side
BP.xposMSmaxfloodlee=xposMSmaxfloodlee; % vertical position max angle flood lee side (after Cisneros et al 2020)
BP.zposMSmaxfloodlee=zposMSmaxfloodlee; % vertical position max angle flood lee side (after Cisneros et al 2020)

BP.xtrough1=xtrough1;                   % x pos trough 1
BP.ytrough1=ytrough1;                   % y pos trough 1
BP.ztrough1=ztrough1;                   % z pos trough 1
BP.xcrest=xcrest;                       % x pos crest
BP.ycrest=ycrest;                       % y pos crest
BP.zcrest=zcrest;                       % depth crest
BP.xtrough2=xtrough2;                   % x pos trough 2
BP.ytrough2=ytrough2;                   % y pos trough 2
BP.ztrough2=ztrough2;                   % z pos trough 2
BP.meanbedformdepth=meanbedformdepth;   % mean bedform depth

BP.MESP=MESP;                           % mean slope ebb steep face
BP.HESP=HESP;                           % mean height ebb steep face
BP.LESP=LESP;                           % mean length ebb steep face
BP.x_esp_beg=x_esp_beg;                 % beginning ebb steep face
BP.x_esp_end=x_esp_end;                 % end ebb steep face
BP.x_esp_begN=x_esp_begN;               % normalised pos beginning ebb steep face
BP.x_esp_endN=x_esp_endN;               % normalised pos end ebb steep face

BP.MFSP=MFSP;                           % mean slope flood steep face
BP.HFSP=HFSP;                           % mean height flood steep face
BP.LFSP=LFSP;                           % mean length flood steep face
BP.x_fsp_beg=x_fsp_beg;                 % beginning flood steep face
BP.x_fsp_end=x_fsp_end;                 % end flood steep face
BP.x_fsp_begN=x_fsp_begN;               % normalised pos beginning flood steep face
BP.x_fsp_endN=x_fsp_endN;               % normalised pos end flood steep face
BP.isthereSP=isthereSP;                 % is there a steep face 0 = no steep face, 1 = ebb steep face, -1 = flood steep face and 2 = both steep face


