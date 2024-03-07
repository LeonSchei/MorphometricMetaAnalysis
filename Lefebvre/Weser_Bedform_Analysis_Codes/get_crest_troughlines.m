function MC=get_crest_troughlines(xr,yr,zr,crestlinethreshold)
%detects bedform crestlines and troughlines
% input data are matrices of xr,yr,zr coordinates with 
% xr longitude/crosswise coordinates
% yr latitude/streamwise coordinates
% zr depth
% 
% The crestlines are detected using a minimum curvature method, vaguely 
% following Ogor (2018)
% The threshold value to detect the crestline can be given. 
% If left blank, a value of the 5th percentile of all minimum curvature
% values is used
%
% output is a big structure with crestline properties, all crestlines (MC.CL)
% and troughlines (MC.TL), and bifurcations (could be improved)
%
% written by Alice Lefebvre, 2021, alefebvre@marum.de

MC=struct; % create output matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% detecting the crests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using the minimum curvature, following Ogor 2018
[~,~,~,P2] = surfature(xr,yr,zr);

if isempty(crestlinethreshold)            % default value if input left blank
    crestlinethreshold=prctile(P2(:),5);       
end
MC.crestlinethreshold=crestlinethreshold;

crestlinesBW=zeros(size(zr));
crestlinesBW(P2<=crestlinethreshold)=1;    % keep only when P2 is less than threshold

% doing some image analysis, feel free to comment/uncomment/play a bit with this
CL=bwmorph(crestlinesBW,'clean');          % Remove isolated pixels (1's surrounded by 0's)
CL=bwmorph(CL,'bridge');                   % Bridge previously unconnected pixels
CL=bwmorph(CL,'fill');                     % Fill isolated interior pixels (0's surrounded by 1's)
CL=bwmorph(CL,'thin',Inf);                 % With N = Inf, remove pixels on the boundaries of objects without allowing objects to break apart

% MC_crtr=double(CL);
% figure
% pcolor(xr,yr,MC_crtr)
% shading flat
% axis equal

% working on bifurcations - that's not very good, it could definitely be improved!
pos_bif_CL_MC=find(bwmorph(CL,'branchpoints')>0);
% for now, only separating bifurcations and saying that each branch is a crestline
for n=1:length(pos_bif_CL_MC)
    CL(pos_bif_CL_MC(n)-2:pos_bif_CL_MC(n))=0;
end
pos_bif2_CL_MC=find(bwmorph(CL,'branchpoints')>0);
for n=1:length(pos_bif2_CL_MC)
    CL(pos_bif2_CL_MC(n)-1:pos_bif2_CL_MC(n))=0;
end

% new analysis now that bifurcations are removed
CL=bwareaopen(CL,8);                       % removing small areas
CLcl = bwconncomp(CL,8);                   % looking for connected points = crestlines

m=0;
% removing those who have the same y but several values in x!!
% probably not very efficient workaround, it could be improved
for ii=1:CLcl.NumObjects
    [~,J]=ind2sub(size(xr),CLcl.PixelIdxList{ii});
    if length(unique(J))<length(J)
        
        pos=CLcl.PixelIdxList{ii};
        [~,IA,~] = unique(J,'stable');
        % keeping only those who have at least 4 points
        if length(pos(IA))>3
            CLcl.PixelIdxList{ii}=pos(IA);     
            m=m+1;
        end
        posn=setdiff(pos,pos(IA));
        % if there is another branch, keeping it if it has more than 3 points
        if isempty(posn)==0
            [~,K]=ind2sub(size(xr),posn);
            if length(unique(K))<length(K)
                [~,IB,~] = unique(K,'stable');
                posn=posn(IB);
            end
            [~,K]=ind2sub(size(xr),posn);
            if length(posn)>3 && mean(diff(K))==1
                CLcl.PixelIdxList{length(CLcl.PixelIdxList)+1}=posn;      % keeping only those who have at least 4 points
                m=m+1;
            end
            clear I J pos posn
        end
    end
end

% updating the matrix
for n=CLcl.NumObjects:-1:1
    if isempty(CLcl.PixelIdxList{n})
        CLcl.PixelIdxList(n)=[];
    end
end
CLcl.NumObjects=length(CLcl.PixelIdxList);

%%%%%%%%%%%%%%%%%%
MC_crtr=zeros(size(zr));
for n=1:CLcl.NumObjects
     MC_crtr(CLcl.PixelIdxList{n})=1;           % BW image of crests
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculating crestline properties
%%%%%% for the orientation, a matrix with the same size as original data and 
% filled with NaNs is created. In places where a crestpoint has been
% calcuated, the orientation is calculated

% crestline orientation
MC.crestO=nan(size(zr));
mid=round(size(xr,2)/2);
for ii=1:CLcl.NumObjects
    xsf=xr(CLcl.PixelIdxList{ii});
    ysf=yr(CLcl.PixelIdxList{ii});
    
    % smoothing the crestline
    ysmooth=smooth(ysf);
    CLcl.smoothed{ii}=ysmooth;

    % finding average waterway point
    [I,~]=ind2sub(size(xr),CLcl.PixelIdxList{ii});
    Li=round(mean(I));
    if Li==length(xr);   Li=Li-1;    end
    
    angledif=nan(length(xsf),1);
    for nn=1:length(xsf)-1
        V2x=xr(Li+1,mid)-xr(Li,mid);
        V2y=yr(Li+1,mid)-yr(Li,mid);
        V1x=xsf(nn+1)-xsf(nn);
        V1y=ysmooth(nn+1)-ysmooth(nn);
        
        angledif(nn)=atan2d(V1x*V2y-V1y*V2x,V1x*V2x+V1y*V2y)-90;
        
    end
    angledif(length(xsf))=angledif(length(xsf)-1);
    angledif(angledif<-90)=angledif(angledif<-90)+180;
    
    MC.crestO(CLcl.PixelIdxList{ii})=angledif;   % putting in a big matrix
    
%     figure
%     line(xsf,ysf,'marker','o','markerfacecolor','k')
%     axis equal
%     line(xsf,ysmooth,'color','b','linewidth',3) 
%     line(xi,yi,'color','b','marker','x') 
%     line(xi(Li),yi(Li),'color','r','marker','o')

end

%%%%%% for the other properties, a matrix with the same size as original data
%  and filled with NaNs is created. 
% the middle point of each crestline is determined, and each crestline
% property is stored there
MC.crest_length=nan(size(zr));
MC.orientation_mean=nan(size(zr));
MC.orientation_std=nan(size(zr));
MC.height_std=nan(size(zr));
MC.NDS=nan(size(zr));              % non-dimensionnal span - if NDS > 1.2, crestline is 3D
MC.perp=0;
% analysis of the crestlines
for n=1:CLcl.NumObjects
    pos=CLcl.PixelIdxList{n};
    POS=pos(round(length(pos)/2));   % middle position
    
    MC.orientation_mean(POS)=nanmean(MC.crestO(pos));
    % if the mean orientation of the crestline is more than 45° compared to
    % the streamwise direction, this is indicated in the matrix MC.prep 
    if find(abs(MC.orientation_mean)>45)
        MC.orientation_mean(POS)=NaN;
        MC.crestO(pos)=NaN;
        MC_crtr(pos)=0;
        MC.perp=MC.perp+1;
    else
        MC.crest_length(POS)=sum(sqrt(diff(xr(pos)).^2+diff(yr(pos)).^2));
        MC.orientation_std(POS)=nanstd(MC.crestO(pos));
        MC.height_std(POS)=nanstd(zr(pos));
    end
    
    Ly=sqrt(diff(xr([pos(1) pos(end)])).^2+diff(yr([pos(1) pos(end)])).^2);
    MC.NDS(POS)=MC.crest_length(POS)/Ly;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% detecting trough points as minimum depth between crest points
TL=zeros(size(zr));
for n=1:size(zr,2)
    d0=MC_crtr(:,n);
    zb=zr(:,n);
    crest_curv_pos=find(d0==1);
    trough_curv_pos=nan(length(crest_curv_pos),1);
    for i=1:length(crest_curv_pos)-1
        trough_curv_pos(i)=crest_curv_pos(i)+find(zb(crest_curv_pos(i):crest_curv_pos(i+1))==min(zb(crest_curv_pos(i):crest_curv_pos(i+1))),1)-1;
    end
    if isempty(crest_curv_pos)==0
        trough_curv_pos(end)=crest_curv_pos(end)+find(zb(crest_curv_pos(end):end)==min(zb(crest_curv_pos(end):end)),1)-1;
    end
    TL(trough_curv_pos,n)=1;
end

% making them into trough lines
TL=bwmorph(TL,'close' );

% doing some image analysis, feel free to comment/uncomment/play a bit with this
TL=bwmorph(TL,'clean');                    % Remove isolated pixels (1's surrounded by 0's)
% TL=bwmorph(TL,'bridge');                 % Bridge previously unconnected pixels
% TL=bwmorph(TL,'fill');                   % Fill isolated interior pixels (0's surrounded by 1's)
TL=bwmorph(TL,'thin',Inf);                 % With N = Inf, remove pixels on the boundaries of objects without allowing objects to break apart

% working on bifurcations - not very good but also not so important for troughs (in the Weser)
pos_bif_TL_MC=find(bwmorph(TL,'branchpoints')>0);       % bifurcations
% for now, only separating bifurcations and saying that each branch is a troughline
for n=1:length(pos_bif_TL_MC)
    TL(pos_bif_TL_MC(n)-2:pos_bif_TL_MC(n))=0;
end
pos_bif2_TL_MC=find(bwmorph(TL,'branchpoints')>0);
for n=1:length(pos_bif2_TL_MC)
    TL(pos_bif2_TL_MC(n)-1:pos_bif2_TL_MC(n))=0;
end

% new analysis now that bifurcations are removed
TL=bwareaopen(TL,4);                       % removing small areas
TLcl = bwconncomp(TL,8);                   % looking for connected points = troughlines

for n=1:TLcl.NumObjects
     MC_crtr(TLcl.PixelIdxList{n})=-1;          % BW image of troughs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% putting in a structure
MC.crtr=MC_crtr;

MC.CL=CLcl;
MC.CL.nb_branch=m;      % just to check how many were with several branches 

MC.TL=TLcl;

MC.pos_bif=nan(size(zr));
MC.pos_bif(pos_bif_CL_MC)=1;

end

% %%%%%%%%%%%%%%%%%%%%% if you want to draw a figure, use lines below
% % % 
% figure
% pcolor(xr,yr,zr)
% axis equal
% shading flat
% % % drawing crestlines
% for n=1:CLcl.NumObjects
%     line(xr(CLcl.PixelIdxList{n}),yr(CLcl.PixelIdxList{n}),'color','k','linewidth',2)
%     %line(xr(CLcl.PixelIdxList{n}),CLcl.smoothed{n},'color','r','linewidth',2)
% end
% % drawing troughlines
% for n=1:TLcl.NumObjects
%     line(xr(TLcl.PixelIdxList{n}),yr(TLcl.PixelIdxList{n}),'color','w','linewidth',2)
% end
% 
% % %drawing crestlines
% % for n=1:CL_ZRcl.NumObjects
% %     line(xr(CL_ZRcl.PixelIdxList{n}),yr(CL_ZRcl.PixelIdxList{n}),'color','k','linewidth',1)
% % end
% % % drawing troughlines
% % for n=1:TL_ZRcl.NumObjects
% %     line(xr(TL_ZRcl.PixelIdxList{n}),yr(TL_ZRcl.PixelIdxList{n}),'color','w','linewidth',2)
% % end
% 
% bifurc=find(MC.pos_bif==1);
% line(xr(bifurc),yr(bifurc),'markerfacecolor','r','marker','o','linestyle','none')
% % 