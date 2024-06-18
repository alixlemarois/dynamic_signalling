function [selected_track_data] = make_biosensor_tracks(im_segmented,nuclei_props,trackfile,selected_track_data)

% Dr Alix LE MAROIS - The Francis Crick Institute - June 2024

% This function uses (x,y,t) coordinates produced by Trackmate to generate
% traces with biosensor and morphological features as computed by the
% 'compute_biosensors_FRET_KTR' function.

% INPUT ARGUMENTS:
% im_segmented : segmented image as produced by the 'compute_biosensors_FRET_KTR' function
% 
% nuclei_props : segmented nuclei properties file as produced by the 'compute_biosensors_FRET_KTR' function
% 
% trackfile : table file from trackmate output. this must be imported into
% matlab and cropped so that the table contains the following variables in
% this order : [track_id quality x_position y_position z_position time time_index].
% note - quality and z_position are not used. 
%
% selected_track_data : if the tracks have already been extracted but the
% biosensor/morphological values have been re-computed, input the old
% 'selected_track_data' file to be overwritten.


% OUTPUT ARGUMENTS:
% selected_track_data : 1-by-nb_of_tracks cell array. Each one contains the
% cellular trace for all biosensor and morphological features.


num_xpix = size(im_segmented,1);
num_ypix = size(im_segmented,2);
if nargin<4
    nb_trackpoints=size(trackfile,1);
    tracknb=1;
    selected_track_data={};
    selected_track_data{1,tracknb}.TRACK_ID(1)=trackfile(1,1);
    selected_track_data{1,tracknb}.xpos(1)=trackfile(1,3)+1;
    selected_track_data{1,tracknb}.ypos(1)=trackfile(1,4)+1;
    selected_track_data{1,tracknb}.t(1)=(trackfile(1,6)+1);
    selected_track_data{1,tracknb}.frame(1)=trackfile(1,7)+1;

    index=1;
    for j=2:nb_trackpoints

        if trackfile(j,1)>trackfile(j-1,1)
            tracknb=tracknb+1;
            index=1;
        else
            index=index+1;
        end

        selected_track_data{1,tracknb}.TRACK_ID(index)=trackfile(j,1);
        selected_track_data{1,tracknb}.xpos(index)=trackfile(j,3)+1;
        selected_track_data{1,tracknb}.ypos(index)=trackfile(j,4)+1;
        selected_track_data{1,tracknb}.t(index)=(trackfile(j,6)+1);
        selected_track_data{1,tracknb}.frame(index)=trackfile(j,7)+1;
    end

    nb_tracks=tracknb;
else
    nb_tracks=max(size(selected_track_data));
    for i = 1:numel(selected_track_data)
        selected_track_data{i}.frame = selected_track_data{i}.t+1;
    end
end

for j=1:nb_tracks

    track_length=max(size(selected_track_data{1,j}.xpos));
    
    for i=1:track_length
              if round(selected_track_data{1,j}.xpos(i))==0
                  selected_track_data{1,j}.xpos(i)=1;
              end
              if round(selected_track_data{1,j}.ypos(i))==0
                  selected_track_data{1,j}.ypos(i)=1;
              end
              if round(selected_track_data{1,j}.xpos(i))>num_xpix
                  selected_track_data{1,j}.xpos(i)=num_xpix;
              end
              if round(selected_track_data{1,j}.ypos(i))>num_ypix
                  selected_track_data{1,j}.ypos(i)=num_ypix;
              end
              [i j];
              round(selected_track_data{1,j}.ypos(i));
              round(selected_track_data{1,j}.xpos(i));
              round(selected_track_data{1,j}.frame(i));
        cellnb=double(im_segmented(round(selected_track_data{1,j}.ypos(i)),round(selected_track_data{1,j}.xpos(i)),round(selected_track_data{1,j}.frame(i))));
    
            if isnan(cellnb);
                cellnb=0;
            end
        
        if cellnb

        selected_track_data{1,j}.ROI_nb(i)=cellnb;
        selected_track_data{1,j}.cell_av_FRET(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.cell_av_FRET(cellnb);
        %selected_track_data{1,j}.KTR_ratio_ring(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.KTR_ratio_ring(cellnb);
        selected_track_data{1,j}.CFP_av(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.CFP_av(cellnb);
        selected_track_data{1,j}.YFP_av(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.YFP_av(cellnb);
        selected_track_data{1,j}.CFP_std(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.CFP_std(cellnb);
        selected_track_data{1,j}.YFP_std(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.YFP_std(cellnb);
        selected_track_data{1,j}.INT_av(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.INT_av(cellnb);
        %selected_track_data{1,j}.KTR_nuc_int(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.KTR_nuc_int(cellnb);
        %selected_track_data{1,j}.KTR_cyto_int(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.KTR_cyto_int(cellnb);
        selected_track_data{1,j}.Area(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.Area(cellnb);
        selected_track_data{1,j}.MajorAxisLength(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.MajorAxisLength(cellnb);
        selected_track_data{1,j}.MinorAxisLength(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.MinorAxisLength(cellnb);
        selected_track_data{1,j}.Eccentricity(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.Eccentricity(cellnb);
        selected_track_data{1,j}.EquivDiameter(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.EquivDiameter(cellnb);
        selected_track_data{1,j}.Perimeter(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.Perimeter(cellnb);
        selected_track_data{1,j}.Aspect_Ratio(i)=nuclei_props{round(selected_track_data{1,j}.frame(i))}.Aspect_Ratio(cellnb);
        
        else
            
        selected_track_data{1,j}.ROI_nb(i)=NaN;
        selected_track_data{1,j}.cell_av_FRET(i)=NaN;
        %selected_track_data{1,j}.KTR_ratio_ring(i)=NaN;
        selected_track_data{1,j}.CFP_av(i)=NaN;
        selected_track_data{1,j}.YFP_av(i)=NaN;
        selected_track_data{1,j}.CFP_std(i)=NaN;
        selected_track_data{1,j}.YFP_std(i)=NaN;
        selected_track_data{1,j}.INT_av(i)=NaN;
        %selected_track_data{1,j}.KTR_nuc_int(i)=NaN;
        %selected_track_data{1,j}.KTR_cyto_int(i)=NaN;
        selected_track_data{1,j}.Area(i)=NaN;
        selected_track_data{1,j}.MajorAxisLength(i)=NaN;
        selected_track_data{1,j}.MinorAxisLength(i)=NaN;
        selected_track_data{1,j}.Eccentricity(i)=NaN;
        selected_track_data{1,j}.EquivDiameter(i)=NaN;
        selected_track_data{1,j}.Perimeter(i)=NaN;
        selected_track_data{1,j}.Aspect_Ratio(i)=NaN;

        end
        
    end
end

end

