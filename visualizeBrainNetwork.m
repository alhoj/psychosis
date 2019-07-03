function h = visualizeBrainNetwork(conMat,thr,clusts,bg,filePrefix)
%conMat - Pre-thresholded connectivity matrix, either weighted or unweighted
%clusts - Spatial maps of nodes in single 3d volume. Each node should have a unique number.
%bg     - Background brain image for visualization
%
% 25.2.2015 Juha Lahnakoski
% V.2.381e-599 pre-alpha

addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/');
addpath('/m/nbe/scratch/braindata/shared/toolboxes/smoothpatch/');
addpath('/m/nbe/scratch/psykoosi/scripts/patchline/')
load('/m/nbe/scratch/psykoosi/masks/brainnetome_EPIgroupMask_inds.mat','regions');

%The x-axis is reversed in this visualization
clusts=clusts(end:-1:1,:,:);
scl=bg.hdr.dime.pixdim(2);
bg=bg.img;
lastClust=max(clusts(:));
conMat(logical(eye(size(conMat))))=0;
find(sum(conMat>thr)>conMat(1,1));
conMatThr=conMat.*(conMat>thr);

% conMat=conMat/max(conMat(:));

for r=1:length(regions.labels)
    figure;
    idx=find(max(conMatThr)>0);
    idx(find(idx<regions.idxStart(r) | idx>regions.idxEnd(r)))=[];

    if ~isempty(idx)
        disp(regions.labels{r})
        
        clusts_temp=clusts;
        for k=1:lastClust
            if ~any(idx==k)
                clusts_temp(clusts_temp==k)=0;
            end
        end

        for k=1:length(idx)
            [x,y,z]=ind2sub(size(clusts_temp),find(clusts_temp==idx(k)));
            centroids(:,idx(k))=mean([x y z]);
        end

%         cmap=jet(lastClust);
        cmap=aaltoColors([1 11 2 9 3 8 6]);
        cmap=max(0,min(1,resample(cmap,lastClust+ceil(lastClust/length(cmap)),length(cmap))));
        cmap(lastClust+1:end,:)=[];
        bgClose=double(imfill(bg>3750,'holes'));

        set(gcf,'color',[1 1 1]);
        set(gca,'visible','off','DataAspectRatioMode','Manual','position',[0 0 1 1]);

        hold on;
        %Smooth and draw background brain
        pp=isosurface(bgClose,.5);
        pp.vertices=pp.vertices/(2/scl);
        pp=smoothpatch(pp,1,4);
        pp=patch(pp,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','faceAlpha',0.3);
        %p=patch(isosurface(bgClose,.5),'edgecolor','none','facecolor',aaltoColors(5)/255,'faceAlpha',.3);
        %isonormals(bgClose,p);

        xlim([-10 size(bg,1)*1.5]);
        ylim([-10 size(bg,2)*1.5]);
        zlim([-10 size(bg,3)*1.5]);
%         xlim([-10 size(bg,1)/(2/scl)]+11);
%         ylim([-10 size(bg,2)/(2/scl)]+11);
%         zlim([-10 size(bg,3)/(2/scl)]+11);
%         xlim([-10 size(bg,1)/2]+11);
%         ylim([-10 size(bg,2)/2]+11);
%         zlim([-10 size(bg,3)/2]+11);

        lighting phong;
        axis equal;

        %Draw ROIs (nodes)
        for k=1:length(idx)
            %p=patch(isosurface(clusts==idx(k),.5),'edgecolor','none','facecolor',cmap(idx(k),:),'faceAlpha',.5);
            %isonormals(clusts==idx(k),p);
            %if any(clusts==idx(k))
            p=isosurface(clusts_temp==idx(k),.5);

            %if ~isempty(p.vertices)
            pSmooth{k}=smoothpatch(p);
            %pClust{k}=
    %         patch(pSmooth{k},'FaceColor',cmap(idx(k),:),'EdgeColor','none','faceAlpha',.35,'AmbientStrength',.5);
            patch(pSmooth{k},'FaceColor',cmap(idx(k),:),'EdgeColor','none','faceAlpha',.5,'AmbientStrength',.5);
            %end;v
            %end;
        end
        %Draw edges
    %     for k=1:length(idx)
    %         for l=k+1:length(idx)
    %             if conMatThr(idx(k),idx(l))>0
    %                 xx=[centroids(2,idx(k)) mean([centroids(2,idx(k)) 45 centroids(2,idx(l))]) centroids(2,idx(l))];
    %                 yy=[centroids(1,idx(k)) mean([centroids(1,idx(k)) 54 centroids(1,idx(l))]) centroids(1,idx(l))];
    %                 zz=[centroids(3,idx(k)) mean([centroids(3,idx(k)) 45 centroids(3,idx(l))]) centroids(3,idx(l))];
    %                 plt=fnplt(cscvn([xx;yy;zz]));
    %                 %pLine{idx(k),idx(l)}=
    %                 patchline(plt(1,:),plt(2,:),plt(3,:),'edgecolor',mean(cmap(idx([k l]),:)),'lineWidth',conMat(idx(k),idx(l))*7);%,'edgealpha',0);%max(cMat(x,y,t),0));
    %                 %plot3([centroids(2,idx(k)) centroids(2,idx(l))],[centroids(1,idx(k)) centroids(1,idx(l))],[centroids(3,idx(k)) centroids(3,idx(l))],'color',mean(cmap(idx([k l]),:)),'linewidth',2);
    %             end
    %         end
    %     end
        %Save different views
        views=[-90 0;180 0;-90 90];
        for curView=1:length(views)
            view(views(curView,1),views(curView,2));
            l=camlight;
            saveas(gcf,sprintf('%s_%s_view_%i.png',filePrefix,regions.labels{r},curView));
            %F=getframe;
            delete(l);
        end
        close all;
        clear p pp
    end
end
end