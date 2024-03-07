function steep_face=tidal_steep_faces(xr,yr,zr,anglethreshold,dirthreshold)
%steep_face=tidal_steep_faces(xr,yr,zr,anglethreshold,dirthreshold)
% creates a structure with steep faces position from 3D bathymetry from a
% tidal environment
% input:
% xr longitude/crosswise coordinates
% yr latitude/streamwise coordinates
% zr depth
% angle threshold to define steep face, usually 15°
% threshold for the direction to define a steep face compared to the main flow direction (usually 45°) 
% 
% ebb is assumed to be towards positive xr
%
% output are given in a structure with 12 matrices (same size as zr) wheer
% position and properties of the steep faces are given
%
% Alice Lefebvre, 2021 alefebvre@marum.de

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position of points where angles are calculated - not sure those are still
% right when using the total angle!
xr_sw=xr(1:end-1,:)+(xr(2:end,:)-xr(1:end-1,:))/2;  xr_sw=[xr_sw;nan(1,size(xr_sw,2))];
yr_sw=yr(1:end-1,:)+(yr(2:end,:)-yr(1:end-1,:))/2;  yr_sw=[yr_sw;nan(1,size(yr_sw,2))];
zr_sw=nan(size(xr_sw));
for i=1:size(zr,2)
    zr_sw(:,i)=interp1(sqrt((xr(:,i)-xr(1,i)).^2+(yr(:,i)-yr(1,i)).^2),zr(:,i),sqrt((xr_sw(:,i)-xr(1,i)).^2+(yr_sw(:,i)-yr(1,i)).^2));
end

% distance between points in the streamwise and crosswise directions
% this is important along a curvilinear grid
swdist=sqrt((xr(2:end,:)-xr(1:end-1,:)).^2+(yr(2:end,:)-yr(1:end-1,:)).^2);
cwdist=sqrt((xr(:,2:end)-xr(:,1:end-1)).^2+(yr(:,2:end)-yr(:,1:end-1)).^2);

% calculating the slopes
gradsw=diff(zr,1,1)./swdist;      gradsw=[gradsw;nan(1,size(gradsw,2))];
gradcw=diff(zr,1,2)./cwdist;      gradcw=[gradcw nan(size(gradcw,1),1)];  

% slope magnitude and direction
[gradd,gradm]=cart2pol(gradcw,gradsw);

slopemag=atand(gradm);                % slope magnitude
slopedir=rad2deg(gradd);              % slope direction, 180 is flood direction; -180 is ebb direction

% figure
% pcolor(xr_sw,yr_sw,slopemag)
% shading flat
% axis equal
% colorbar

% removing those which are not ebb or flood direction, using the direction threshold
slopemag(slopedir<-90-dirthreshold)=NaN;
slopemag(slopedir>90+dirthreshold)=NaN;
slopemag(slopedir>-90+dirthreshold & slopedir<90-dirthreshold)=NaN;

% putting ebb as negative slopes
slopemag(atand(gradsw)<0)=-slopemag(atand(gradsw)<0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% finding steep faces
slopetoanalyse=slopemag;    % here also possibility to do the analysis on another slope, e.g. streamwise=atand(gradsw)

% ebb steep face
ebbslopeBW=zeros(size(slopetoanalyse));
ebbslopeBW(slopetoanalyse<=-anglethreshold)=1;          % keep only slope steeper than 15°
BW_swebb=im2double(ebbslopeBW,'indexed');               % BW image of crests
BW2_swebb = bwareaopen(BW_swebb,4);                     % removing small areas
steep_face_sw_ebb=im2double(BW2_swebb,'indexed');
steep_face_sw_ebb(steep_face_sw_ebb==0)=-1;
steep_face_sw_filled_ebb=steep_face_sw_ebb;
effect_ebb=0;

% here doing a steep_face_sw_filled where the "holes" are filled, that's
% old, there is probably something better to be done
for n=1:size(steep_face_sw_ebb,1)
    MAT = findseq(steep_face_sw_filled_ebb(n,:));
    MAT(MAT(:,1)==-1,:)=[];
    
    % remove holes that are less than 4 cells
    dm4=find((MAT(2:end,2)-MAT(1:end-1,3))<6);
    effect_ebb=effect_ebb+length(dm4);
    for i=1:length(dm4)
        steep_face_sw_filled_ebb(n,MAT(dm4(i),3):MAT(dm4(i)+1,2))=1;
    end
end

steep_face_sw_ebb=steep_face_sw_filled_ebb;    % here I decided to keep the filled steep face for the analysis

% flood steep face
floodslopeBW=zeros(size(slopetoanalyse));
floodslopeBW(slopetoanalyse>=anglethreshold)=1;             % keep only slope steeper than 15°
BW_swflood=im2double(floodslopeBW,'indexed');               % BW image of crests
BW2_swflood = bwareaopen(BW_swflood,4);                      % removing small areas
steep_face_sw_flood=im2double(BW2_swflood,'indexed');
steep_face_sw_flood(steep_face_sw_flood==0)=-1;
steep_face_sw_filled_flood=steep_face_sw_flood;
effect_flood=0;
% here doing a steep_face_sw_filled where the "holes" are filled
for n=1:size(steep_face_sw_ebb,1)
    MAT = findseq(steep_face_sw_filled_flood(n,:));
    MAT(MAT(:,1)==-1,:)=[];
    
    % remove holes that are less than 4 cells
    dm4=find((MAT(2:end,2)-MAT(1:end-1,3))<6);
    effect_flood=effect_flood+length(dm4);
    for i=1:length(dm4)
        steep_face_sw_filled_flood(n,MAT(dm4(i),3):MAT(dm4(i)+1,2))=1;
    end
end

steep_face_sw_flood=steep_face_sw_filled_flood;    % here I decided to keep the filled steep face for the analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% surface area
% surface area of each steep face, quite simple, I guess this could be improved
% ebb 
SA_ebb=nan(size(zr));
esf=find(steep_face_sw_ebb>0);
SA_ebb(esf)=swdist(esf)./cosd(slopetoanalyse(esf));

% flood
SA_flood=nan(size(zr));
fsf=find(steep_face_sw_flood>0);
SA_flood(fsf)=swdist(fsf)./cosd(slopetoanalyse(fsf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% orientation
% calculate the steep face orientation - ebb
steep_face.ebb_O=nan(size(zr));
spe=zeros(size(zr_sw));  spe(steep_face_sw_ebb==1)=1;
CC_ebb=bwconncomp(spe,8);
mid=round(size(xr,2)/2);
for ii=1:CC_ebb.NumObjects

    xsf=xr_sw(CC_ebb.PixelIdxList{ii});
    ysf=yr_sw(CC_ebb.PixelIdxList{ii});
    
    % doing a smoothingsline regression through the points
    [curve]=fit(xsf,ysf,'smoothingspline','SmoothingParam',0.005);
    xsm=unique(ceil(xsf));
    ysm=curve(xsm);

    % finding average waterway point
    [I,~]=ind2sub(size(xr_sw),CC_ebb.PixelIdxList{ii});
    Li=round(mean(I));
    if Li==length(xr_sw);   Li=Li-1;    end
    
    angledif=nan(length(xsm)-1,1);
    posdif=nan(length(xsm)-1,1);
    
    for nn=1:length(xsm)-1
        V2x=xr_sw(Li+1,mid)-xr_sw(Li,mid);
        V2y=yr_sw(Li+1,mid)-yr_sw(Li,mid);
        V1x=xsm(nn+1)-xsm(nn);
        V1y=ysm(nn+1)-ysm(nn);
        
        angledif(nn)=atan2d(V1x*V2y-V1y*V2x,V1x*V2x+V1y*V2y)-90;
        
        n1=find(abs(yr-ysm(nn))==min(abs(yr-ysm(nn))));
        %n2=find(abs(xr(n1)-xsm(nn))==min(abs(xr(n1)-xsm(nn))));
        posdif(nn)=n1(find(abs(xr(n1)-xsm(nn))==min(abs(xr(n1)-xsm(nn))),1));
    end
    angledif(angledif<-90)=angledif(angledif<-90)+180;
    
    steep_face.ebb_O(posdif)=angledif;   % putting in a big matrix
    
   
%     figure
%     line(xsf,ysf,'marker','o','markerfacecolor','k')
%     axis equal
%     line(xsm,ysm,'color','b','linewidth',3) 
%     line(xr(:,mid),yr(:,mid),'color','c','marker','.')
%     line(xr(round(mean(I)),mid),yr(round(mean(I)),mid),'color','k','marker','o')
%     line(xr_sw(posdif),yr_sw(posdif),'color','r','marker','.')
end
%
% calculate the steep face orientation - flood
steep_face.flood_O=nan(size(zr));
spf=zeros(size(zr_sw));  spf(steep_face_sw_flood==1)=1;
CC_flood=bwconncomp(spf,8);

for ii=1:CC_flood.NumObjects
    xsf=xr_sw(CC_flood.PixelIdxList{ii});
    ysf=yr_sw(CC_flood.PixelIdxList{ii});
    
    % doing a smoothingsline regression through the points
    [curve]=fit(xsf,ysf,'smoothingspline','SmoothingParam',0.005);
    xsm=unique(ceil(xsf));
    ysm=curve(xsm);

    % finding average waterway point
    [I,~]=ind2sub(size(xr_sw),CC_flood.PixelIdxList{ii});
    Li=round(mean(I));
    if Li==length(xr_sw);   Li=Li-1;    end
    
    angledif=nan(length(xsm)-1,1);
    posdif=nan(length(xsm)-1,1);
    
    for nn=1:length(xsm)-1
        V2x=xr_sw(Li+1,mid)-xr_sw(Li,mid);
        V2y=yr_sw(Li+1,mid)-yr_sw(Li,mid);
        V1x=xsm(nn+1)-xsm(nn);
        V1y=ysm(nn+1)-ysm(nn);
        
        angledif(nn)=atan2d(V1x*V2y-V1y*V2x,V1x*V2x+V1y*V2y)-90;
        
        n1=find(abs(yr-ysm(nn))==min(abs(yr-ysm(nn))));
        %n2=find(abs(xr(n1)-xsm(nn))==min(abs(xr(n1)-xsm(nn))));
        posdif(nn)=n1(find(abs(xr(n1)-xsm(nn))==min(abs(xr(n1)-xsm(nn))),1));
    end
    angledif(angledif<-90)=angledif(angledif<-90)+180;
    
    steep_face.flood_O(posdif)=angledif;   % putting in a big matrix
    
%     figure
%     line(xsf,ysf,'marker','o','markerfacecolor','k')
%     axis equal
%     line(xsm,ysm,'color','b','linewidth',3) 
%     line(xr(:,mid),yr(:,mid),'color','c','marker','.')
%     line(xr(round(mean(I)),mid),yr(round(mean(I)),mid),'color','k','marker','o')
%     line(xr(posdif),yr(posdif),'color','r','marker','.')
end

%%%% putting results in a matrix, -1 are flood steep faces, 1 are ebb steep faces
allsteepface=zeros(size(zr_sw));
allsteepface(steep_face_sw_flood==1)=-1;
allsteepface(steep_face_sw_ebb==1)=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slopemag=atand(slopemag);     % again, just all the magnitude, without those which were removed
slopemag(atand(gradsw)<0)=-slopemag(atand(gradsw)<0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output structure   
steep_face.slopemag=slopemag;                % slope magnitude
steep_face.slopedir=slopedir;                % slope direction
steep_face.xr_sw=xr_sw;                      % x position of angles
steep_face.yr_sw=yr_sw;                      % y position of angles
steep_face.zr_sw=zr_sw;                      % z value at angle position
steep_face.ebb=steep_face_sw_ebb;            % positions of ebb steep faces
steep_face.flood=steep_face_sw_flood;        % position of flood steep faces
steep_face.both=allsteepface;                % position of all steep faces (-1 flood, 1 ebb)
steep_face.SA_ebb=SA_ebb;                    % ebb surface area
steep_face.SA_flood=SA_flood;                % flood surface area