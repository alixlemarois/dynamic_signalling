function [nuclei_props, cell_av_FRET_image, t] = ...
compute_biosensors(im_path, isPh, isKTR, area_min, thresholds, segmentation_type, core_size, data_order)

% Dr Alix LE MAROIS - The Francis Crick Institute - June 2024

% This function uses raw CFP and YFP FRET images as well as KTR images
% along with matching nuclear segmentation masks to generate false-colour
% FRET images as well as tables containing FRET / KTR biosensor morphological values for each
% cell and time point. 

% INPUT ARGUMENTS:
% im_path : path to the directory where the raw and segmented images are
% found. The raw images MUST have 'raw.tif' at the end of their file name,
% and the segmented images MUST match the raw image name but be terminated
% by 'segmented.tif' instead of 'raw.tif'.
% isPh : does the image contain a phase channel. 0 if not, 1 if yes.
% isKTR : does the image contain a KTR sensor channel. 0 if not, 1 if yes.
% area_min : threshold area for a segmented spot to be considered a cell
% (in pixels).
% thresholds : minimum intensity thresholds for both CFP+YFP and KTR images
% - 2-by-1 matrix.
% segmentation_type : segmentation algorithm used in fiji to generate
% segmentation masks. 'ilastik' or 'stardist'.
% core_size : size of nuclear core square area for CN ratio quantification.
% data_order : order in which channels and time are indexed in the file. 
% If indexing over time first, 'T_C' (ie CFP_1 ... CFP_n ; YFP_1 ... YFP_n ; mRuby_1 ... mRuby_n).
% If indexing over channels first, 'C_T' (ie CFP1, YFP_1, mRuby_1 ; ... ;
% CFPn, YFP_n, mRuby_n).

% OUTPUT ARGUMENTS:
% All output files are saved by default in the working directory. If
% needed, variables for individual files can be written into the Matlab
% workspace. However, these will be overwritten each time a new file from
% the directory is imported and analysed.

files=dir(im_path)
intmax=4095;

for i=1:size(files,1)   
    i
    if strfind(files(i).name, '_raw.tif');                
        image=bfopen(files(i).name);
            %% populate raw image cubes with data
            series=image{1,1};

            nb_xpix=size(series{1,1},1)
            nb_ypix=size(series{1,1},2)

            if isPh==0
                if isKTR==0;
                    nb_timeframes=length(series)/2
                elseif isKTR==1;
                    nb_timeframes=length(series)/3
                end
            elseif isPh==1;
                if isKTR==0;
                    nb_timeframes=length(series)/3
                elseif isKTR==1;
                    nb_timeframes=length(series)/4
                end
            end            
            
            CFP_raw=zeros(nb_xpix,nb_ypix,nb_timeframes);
            YFP_raw=zeros(nb_xpix,nb_ypix,nb_timeframes);
            KTR_raw=zeros(nb_xpix,nb_ypix,nb_timeframes);
            phase_im=zeros(nb_xpix,nb_ypix,nb_timeframes);
            CFP_mask_sat=zeros(nb_xpix,nb_ypix,nb_timeframes);
            YFP_mask_sat=zeros(nb_xpix,nb_ypix,nb_timeframes);
            KTR_mask_sat=zeros(nb_xpix,nb_ypix,nb_timeframes);
            
           
            for t=1:nb_timeframes
                if isPh==0;
                    if isKTR==0;
                        switch data_order
                            case 'T_C'
                                CFP_raw(:,:,t)=double(series{t,1});
                                YFP_raw(:,:,t)=double(series{nb_timeframes+t,1});
                                CFP_mask_sat(:,:,t)=CFP_raw(:,:,t)<intmax;
                                YFP_mask_sat(:,:,t)=YFP_raw(:,:,t)<intmax;
                            case 'C_T'
                                CFP_raw(:,:,t)=double(series{2*t-1,1});
                                YFP_raw(:,:,t)=double(series{2*t,1});
                                CFP_mask_sat(:,:,t)=CFP_raw(:,:,t)<intmax;
                                YFP_mask_sat(:,:,t)=YFP_raw(:,:,t)<intmax;
                        end
                    elseif isKTR==1;
                        switch data_order
                            case 'T_C' % file contains all CFP, then all YFP, then all mRuby - when FV3000 images exported in FV software
                                YFP_raw(:,:,t)=double(series{nb_timeframes+t,1});
                                KTR_raw(:,:,t)=double(series{2*nb_timeframes+t,1});
                                CFP_raw(:,:,t)=double(series{t,1});
                                CFP_mask_sat(:,:,t)=CFP_raw(:,:,t)<intmax;
                                YFP_mask_sat(:,:,t)=YFP_raw(:,:,t)<intmax;
                                KTR_mask_sat(:,:,t)=KTR_raw(:,:,t)<intmax;
                            case 'C_T' % file contains all channels at t1 then t2, etc. - when images are generated from raw via imageJ
                                YFP_raw(:,:,t)=double(series{3*t-1,1});
                                KTR_raw(:,:,t)=double(series{3*t,1});
                                CFP_raw(:,:,t)=double(series{3*t-2,1});
                                CFP_mask_sat(:,:,t)=CFP_raw(:,:,t)<intmax;
                                YFP_mask_sat(:,:,t)=YFP_raw(:,:,t)<intmax;
                                KTR_mask_sat(:,:,t)=KTR_raw(:,:,t)<intmax;
                        end
                        t
                    end
                elseif isPh==1;
                   if isKTR==0;
                        CFP_raw(:,:,t)=double(series{3*t-2,1});
                        YFP_raw(:,:,t)=double(series{3*t-1,1});
                        phase_im(:,:,t)=double(series{3*t,1});
                        CFP_mask_sat(:,:,t)=CFP_raw(:,:,t)<intmax;
                        YFP_mask_sat(:,:,t)=YFP_raw(:,:,t)<intmax;
                    elseif isKTR==1;
                        CFP_raw(:,:,t)=double(series{4*t-3,1});
                        YFP_raw(:,:,t)=double(series{4*t-2,1});           
                        KTR_raw(:,:,t)=double(series{4*t-1,1});
                        phase_im(:,:,t)=double(series{4*t,1})*10;
                        CFP_mask_sat(:,:,t)=CFP_raw(:,:,t)<intmax;
                        YFP_mask_sat(:,:,t)=YFP_raw(:,:,t)<intmax;
                        KTR_mask_sat(:,:,t)=KTR_raw(:,:,t)<intmax;
                    end

                end
            end
            
            
            %% subtract background

            for k=1:size(files,1)
                %if strcmp(strcat(extractBetween(files(i).name, 1,".oir"),"_ilastik_segmentation.mat"), files(k).name)
                if strcmp(strcat(extractBetween(files(i).name, 1,"_raw.tif"),"_segmented.tif"), files(k).name)                       
                    segm_output=bfopen(files(k).name);      
                    segm_output=segm_output{1,1};
                end
            end
                
            CFP_bgc=zeros(nb_xpix,nb_ypix,nb_timeframes);
            YFP_bgc=zeros(nb_xpix,nb_ypix,nb_timeframes);
            KTR_bgc=zeros(nb_xpix,nb_ypix,nb_timeframes);
            bg_CFP=zeros(1,nb_timeframes);
            bg_YFP=zeros(1,nb_timeframes);
            bg_KTR=zeros(1,nb_timeframes);
        
            cell_av_FRET_image=zeros(nb_xpix,nb_ypix,nb_timeframes);
            im_segmented=zeros(nb_xpix,nb_ypix,nb_timeframes);
            nuclei_props={}; 
            pop_statistics=table();
            
            se=strel('disk',20,0);
            
            for t=1:nb_timeframes
            %for t = 1:5
                t
                switch segmentation_type
                    case 'stardist'
                        segmented = double(segm_output{t,1});
                        labels = unique(segmented);
                        labels(labels == 0)=[];
                        new_label = 1
                        for l = 1:numel(labels);
                            label_mask = segmented == labels(l) ;
                            im_segmented(:,:,t) = im_segmented(:,:,t) + label_mask.*new_label;
                            new_label = new_label+1;
                        end
                        bg_mask=~imdilate(im_segmented(:,:,t)>0,se);

                        
                    case 'ilastik'
                        cell_code = 2 ;                      
                        cells=(segm_output{t,1}==cell_code);
                        cellfill=imfill(cells, 'holes');
                        cell_nosmall=bwareaopen(cellfill, 10,4);  
                        im_segmented(:,:,t)=cell_nosmall;
                        im_segmented(:,:,t)=bwlabel(im_segmented(:,:,t),4);
                        bg_mask=imdilate(im_segmented(:,:,t)==1,se);
                        
                    case 'cellpose'
                        cells=segm_output{t,1};
                        im_segmented(:,:,t)=cells;
                        bg_mask=~imdilate(im_segmented(:,:,t)>0,se);
                        
                end              
            
                CFP_bg=reshape(CFP_raw(:,:,t).*bg_mask,nb_xpix*nb_ypix,1);
                YFP_bg=reshape(YFP_raw(:,:,t).*bg_mask,nb_xpix*nb_ypix,1);
                
                CFP_bg(CFP_bg==0)=[];
                YFP_bg(YFP_bg==0)=[];
            
                bg_CFP(t)=mean(CFP_bg,'omitnan');
                bg_YFP(t)=mean(YFP_bg,'omitnan');
            
                CFP_bgc(:,:,t)=CFP_raw(:,:,t)-bg_CFP(t);
                YFP_bgc(:,:,t)=YFP_raw(:,:,t)-bg_YFP(t);
                
                % NB no background subtraction for KTR image for now.
                KTR_bgc(:,:,t) = KTR_raw(:,:,t);

                %% generate table of nuclear properties
                warning('off', 'all');
                nuclei_props{t}=regionprops('table',im_segmented(:,:,t),'Area','Centroid', 'Eccentricity','EquivDiameter' ,'MajorAxisLength','MinorAxisLength','Perimeter');
                nuclei_props{t}(nuclei_props{t}.Area==0,:)=[];
                nuclei_props{t}=table2array(nuclei_props{t});
                nuclei_props{t}=array2table(nuclei_props{t});

                if height(nuclei_props{t})>0
                    nuclei_props{t}.Properties.VariableNames={'Area','Centroid_x','Centroid_y','MajorAxisLength','MinorAxisLength','Eccentricity','EquivDiameter','Perimeter'};
                    nuclei_props{t}.Aspect_Ratio=nuclei_props{t}.MajorAxisLength./nuclei_props{t}.MinorAxisLength;       
                    nuclei_props{t}.ROI_number=linspace(1,height(nuclei_props{t}),height(nuclei_props{t}))';

                else,
                    nuclei_props{t}=NaN;
                end

                %% FRET RATIO computation
                if size(nuclei_props{t})>=1
                    % find CFP and YFP pixels - above background and
                    % non-saturated ### ONLY PIXELS fulfilling conditions
                    % in BOTH CHANNELS ###
                    
                    CFP_props=regionprops('table',im_segmented(:,:,t), YFP_mask_sat(:,:,t).*CFP_mask_sat(:,:,t).*CFP_bgc(:,:,t).*(CFP_bgc(:,:,t)>thresholds(1)).*(YFP_bgc(:,:,t)>thresholds(1)),'PixelValues');
                    CFP_props.Properties.VariableNames={'CFP_pixelvalues'};
                    YFP_props=regionprops('table',im_segmented(:,:,t), CFP_mask_sat(:,:,t).*YFP_mask_sat(:,:,t).*YFP_bgc(:,:,t).*(YFP_bgc(:,:,t)>thresholds(1)).*(CFP_bgc(:,:,t)>thresholds(1)),'PixelValues');
                    YFP_props.Properties.VariableNames={'YFP_pixelvalues'};
                    

                    FRET_props=[CFP_props YFP_props];

                    for k=1:height(FRET_props)
                        FRET_props.CFP_pixelvalues{k,1}(FRET_props.CFP_pixelvalues{k,1}==0)=[];
                        FRET_props.YFP_pixelvalues{k,1}(FRET_props.YFP_pixelvalues{k,1}==0)=[];
                            if (~isempty(FRET_props.CFP_pixelvalues{k,1}) && ~isempty(FRET_props.YFP_pixelvalues{k,1}))
                                 nuclei_props{t}.CFP_av(k)=mean(FRET_props.CFP_pixelvalues{k,1},'omitnan');
                                 nuclei_props{t}.YFP_av(k)=mean(FRET_props.YFP_pixelvalues{k,1},'omitnan');
                                 nuclei_props{t}.CFP_std(k)=std(FRET_props.CFP_pixelvalues{k,1},'omitnan');
                                 nuclei_props{t}.YFP_std(k)=std(FRET_props.YFP_pixelvalues{k,1},'omitnan');
                            else,
                                 nuclei_props{t}.CFP_av(k)=NaN;
                                 nuclei_props{t}.YFP_av(k)=NaN;
                                 nuclei_props{t}.CFP_std(k)=NaN;
                                 nuclei_props{t}.YFP_std(k)=NaN;
                            end
                        end
                end
                
                %%
                
                if isKTR
                    
                    if numel(nuclei_props{t})>1
                    
                        % compute KTR ratio

                        % create square around nuclear centroids
                        % square size = 2*core_size+1.

                        centroid_mask=zeros(nb_xpix, nb_ypix);             

                        for m=1:height(nuclei_props{t})
                            centroid_mask(...
                                max([round(nuclei_props{t}.Centroid_y(m))-core_size,1]):min([round(nuclei_props{t}.Centroid_y(m))+core_size nb_ypix]),...
                                max([round(nuclei_props{t}.Centroid_x(m))-core_size,1]):min([round(nuclei_props{t}.Centroid_x(m))+core_size nb_xpix])...
                                )=1;
                        end

                        %dilate nuclei and create ring mask
                        se_ring=strel('disk',2, 8);
                        nuclei_dil_mask=imdilate(im_segmented(:,:,t),se_ring);
                        nuclei_ring_mask=nuclei_dil_mask-im_segmented(:,:,t);

                        % display module for debugging
        %                 figure
        %                 subplot(1,3,1)
        %                 imshow(KTR_bgc(:,:,t).*centroid_mask)
        %                 caxis auto
        %                 subplot(1,3,2)
        %                 imshow(KTR_bgc(:,:,t).*(nuclei_ring_mask>0))
        %                 caxis auto
        %                 subplot(1,3,3)
        %                 imshow(KTR_bgc(:,:,t))
        %                 caxis auto
        %                 aaaa

                        KTR_cyto_props=regionprops('table',nuclei_ring_mask, KTR_bgc(:,:,t).*(KTR_bgc(:,:,t)>0).*KTR_mask_sat(:,:,t),'PixelValues');
                        KTR_cyto_props.Properties.VariableNames={'KTR_cyto_int'};

                        for p=1:height(KTR_cyto_props)
                                % remove zero valued pixels in cyto ring mask
                                KTR_cyto_props.KTR_cyto_int{p,1}(KTR_cyto_props.KTR_cyto_int{p,1}==0)=[];
                             if ~isnan(KTR_cyto_props.KTR_cyto_int{p,1})
                                % compute mean intensity of upper 50% cyto ring pixels.
                                nuclei_props{t}.KTR_cyto_int(p)=mean(KTR_cyto_props.KTR_cyto_int{p,1}(KTR_cyto_props.KTR_cyto_int{p,1}>=median(KTR_cyto_props.KTR_cyto_int{p,1})),'omitnan');
                                %nuclei_props.KTR_cyto_int(i)=mean(cyto_props.KTR_cyto_int{i,1},'omitnan');
                             else,
                                nuclei_props{t}.KTR_cyto_int(p)=NaN;
                             end
                        end

                        centroid_mask=centroid_mask(1:nb_xpix,1:nb_ypix);
                        KTR_nuc_props=regionprops('table',im_segmented(:,:,t),KTR_bgc(:,:,t).*centroid_mask.*(KTR_bgc(:,:,t)>0).*KTR_mask_sat(:,:,t),'PixelValues');
                        KTR_nuc_props.Properties.VariableNames={'KTR_nuc_int'};

                        for q=1:height(KTR_nuc_props)
                                % remove zero valued pixels in nuclear mask
                                KTR_nuc_props.KTR_nuc_int{q,1}(KTR_nuc_props.KTR_nuc_int{q,1}==0)=[];
                             if ~isnan(KTR_nuc_props.KTR_nuc_int{q,1})
                                 % compute mean nuclear intensity value.
                                nuclei_props{t}.KTR_nuc_int(q)=mean(KTR_nuc_props.KTR_nuc_int{q,1},'omitnan');
                             else,
                                nuclei_props{t}.KTR_nuc_int(q)=NaN;
                             end
                        end

                        %compute KTR ratio
                        nuclei_props{t}.KTR_ratio_ring=nuclei_props{t}.KTR_cyto_int./(nuclei_props{t}.KTR_cyto_int+nuclei_props{t}.KTR_nuc_int);

                        %remove sub-threshold and sub-area nuclei
                        nuclei_props{t}.KTR_ratio_ring(nuclei_props{t}.Area<area_min)=NaN;
                        nuclei_props{t}.KTR_ratio_ring(nuclei_props{t}.KTR_cyto_int+nuclei_props{t}.KTR_nuc_int<thresholds(2))=NaN;
                    
                    end
                end
                          
                
                if find(im_segmented(:,:,t));
                    
                    nuclei_props{t}.cell_av_FRET=nuclei_props{t}.YFP_av./nuclei_props{t}.CFP_av; 
                    nuclei_props{t}.INT_av=nuclei_props{t}.YFP_av+nuclei_props{t}.CFP_av; 
                    nuclei_props{t}.cell_av_FRET(nuclei_props{t}.Area<area_min)=NaN;
                    nuclei_props{t}.cell_av_FRET(nuclei_props{t}.CFP_av<thresholds(1))=NaN;
                    nuclei_props{t}.cell_av_FRET(nuclei_props{t}.YFP_av<thresholds(1))=NaN;
                    
                    if isKTR ==1
                        nuclei_props{t}.KTR_ratio_ring(nuclei_props{t}.KTR_cyto_int+nuclei_props{t}.KTR_nuc_int<thresholds(2))=NaN;
                    end
                    
                    nuclei_props{t}.KTR_ratio_ring(nuclei_props{t}.Area<area_min)=NaN;
                    %nuclei_props{t}=[nuclei_props{t} FRET_props];

                    % compute FRET image
                    for n=1:height(nuclei_props{t});
                        %nuclei_props{t}.ROI_number(n)=im_segmented(round(nuclei_props{t}.Centroid_y(n)),round(nuclei_props{t}.Centroid_x(n)),t);
                        mask_segm=(im_segmented(:,:,t)==nuclei_props{t}.ROI_number(n));
                        if ~isnan(nuclei_props{t}.cell_av_FRET(n));
                            cell_av_FRET_image(:,:,t)=cell_av_FRET_image(:,:,t)+(mask_segm*nuclei_props{t}.cell_av_FRET(n));
                        end
                    end
                
                    pop_statistics.mean_FRET(t)=mean(nuclei_props{t}.cell_av_FRET,'omitnan');
                    pop_statistics.std_FRET(t)=std(nuclei_props{t}.cell_av_FRET,'omitnan');
                    pop_statistics.mean_YFP(t)=mean(nuclei_props{t}.YFP_av,'omitnan');
                    pop_statistics.mean_CFP(t)=mean(nuclei_props{t}.CFP_av,'omitnan');
                    pop_statistics.std_YFP(t)=std(nuclei_props{t}.YFP_av,'omitnan');
                    pop_statistics.std_CFP(t)=std(nuclei_props{t}.CFP_av,'omitnan');
                    pop_statistics.cellnumber(t)=height(nuclei_props{t});
                    pop_statistics.bg_CFP(t)=bg_CFP(t);
                    pop_statistics.bg_YFP(t)=bg_YFP(t); 
                    pop_statistics.time_frame(t)=t;
                    
                    if isKTR == 1;
                        pop_statistics.bg_KTR(t)=bg_KTR(t);
                        pop_statistics.KTR_nuc_int_mean(t)=mean(nuclei_props{t}.KTR_nuc_int,'omitnan');
                        pop_statistics.KTR_cyto_int_mean(t)=mean(nuclei_props{t}.KTR_cyto_int,'omitnan');
                        pop_statistics.KTR_ratio_mean(t)=mean(nuclei_props{t}.KTR_ratio_ring,'omitnan');
                        pop_statistics.KTR_ratio_std(t)=std(nuclei_props{t}.KTR_ratio_ring,'omitnan');
                    end
                    
                end
                
                
    
            end
        save(strcat(extractBetween(files(i).name,1,"_raw.tif"),"_segmented_props.mat"),'nuclei_props', 'cell_av_FRET_image', 'im_segmented','pop_statistics');
        if isKTR && isPh
            save(strcat(extractBetween(files(i).name,1,"_raw.tif"),"_preprocess.mat"),'CFP_bgc','YFP_bgc','phase_im','KTR_bgc');
        elseif isPh
            save(strcat(extractBetween(files(i).name,1,"_raw.tif"),"_preprocess.mat"),'CFP_bgc','YFP_bgc','phase_im');
        elseif isKTR
            save(strcat(extractBetween(files(i).name,1,"_raw.tif"),"_preprocess.mat"),'CFP_bgc','YFP_bgc','KTR_bgc');
        else,
           save(strcat(extractBetween(files(i).name,1,"_raw.tif"),"_preprocess.mat"),'CFP_bgc','YFP_bgc');
        end
end
end