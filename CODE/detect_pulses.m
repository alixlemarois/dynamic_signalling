function [track_metrics, track_metrics_summary_table] = detect_pulses(trackfile, exclude_mitosis, frame_rate, data_type, smoothing_factor, min_peak_prominence)
    
if nargin<6
    min_peak_prominence = 0.05
end

if nargin<5
    smoothing_factor = 4;
end

if nargin<4
    data_type = 'trackmate';
end
  
nb_tracks = max(size(trackfile))
track_metrics_summary_table = table();
        
for i = 1:nb_tracks

    % remove mitotic segments of track by replacing with NaN's.
    if exclude_mitosis            
        i                
         % case 1 : manual tracking or trackmate with mitotic track
         if strcmp(data_type, 'manual') || (strcmp(data_type, 'trackmate') && trackfile{i}.ismitotic)
             % if tracks not previously aligned
             if ~isfield(trackfile{i},'t_aligned')                     
                 [trackfile{i}.t_aligned]=align_tracks(trackfile{i});
             end
             % remove data points 2h before mitosis
             mitotic_start = max([find(trackfile{i}.t_aligned == round(-2*60/frame_rate)) 1]);
             mitotic_end = find(trackfile{i}.t_aligned == 0);
             trackfile{i}.cell_av_FRET_nomitosis = trackfile{i}.cell_av_FRET(:,1) ;
             trackfile{i}.cell_av_FRET_nomitosis(mitotic_start:mitotic_end) = NaN ;   
                          % case 2 : trackmate tracking with mitotic tracks
         elseif strcmp(data_type, 'trackmate_xml') && ismember(trackfile{i}.fate, [2 3]);
             if isempty(trackfile{i}.mitosis_time_abs)                   
                 [trackfile{i}.t_aligned]=align_tracks(trackfile{i}.track_data);
             else,                     
                 trackfile{i}.t_aligned = trackfile{i}.track_data.tpos'-trackfile{i}.mitosis_time_abs(1);
             end

              % remove data points 2h before mitosis
              mitotic_start = max([find(trackfile{i}.t_aligned == round(-2*60/frame_rate)) 1]);
              mitotic_end = find(trackfile{i}.t_aligned == 0);
              trackfile{i}.cell_av_FRET_nomitosis = trackfile{i}.cell_av_FRET(:,1) ;
              trackfile{i}.cell_av_FRET_nomitosis(mitotic_start:mitotic_end) = NaN ; 

         else
             trackfile{i}.cell_av_FRET_nomitosis = trackfile{i}.cell_av_FRET(:,1) ;
             trackfile{i}.t_aligned = trackfile{i}.track_data.tpos';
         end                                                  
    else
        trackfile{i}.cell_av_FRET_nomitosis = trackfile{i}.cell_av_FRET;
        trackfile{i}.t_aligned = trackfile{i}.frame;
    end                       

    % compute moving average
    trackfile{i}.cell_av_FRET_smooth = movmean(trackfile{i}.cell_av_FRET_nomitosis, smoothing_factor, 'omitnan');
    track_metrics{i}.cell_av_FRET_smooth = trackfile{i}.cell_av_FRET_smooth;
    track_metrics{i}.t_aligned = trackfile{i}.t_aligned;
    track_metrics{i}.track_nb = i ;
    track_metrics{i}.peak_metrics=table();

    switch data_type        
        case 'trackmate'
            time_vector = trackfile{i}.frame*frame_rate/60;                  
        case 'manual'
            time_vector = trackfile{i}.t*frame_rate/60;
        case 'trackmate_xml'
            time_vector = trackfile{i}.track_data.tpos'*frame_rate/60;
    end

    if max(size(time_vector))>=3
         % detect peaks with minimum prominence of minpeakprominence    
        [pks,locs,widths,proms] = findpeaks(trackfile{i}.cell_av_FRET_smooth, time_vector, 'MinPeakProminence', min_peak_prominence);

        track_metrics{i}.peak_metrics.peak_location=locs';
        track_metrics{i}.peak_metrics.peak_height=pks';
        track_metrics{i}.peak_metrics.peak_width=widths';
        track_metrics{i}.peak_metrics.peak_prominence=proms';
        track_metrics{i}.peak_metrics.peak_baseline=track_metrics{i}.peak_metrics.peak_height-track_metrics{i}.peak_metrics.peak_prominence;

        track_metrics{i}.peak_frequency = height(track_metrics{i}.peak_metrics)./((numel(trackfile{i}.cell_av_FRET_smooth)-sum(isnan(trackfile{i}.cell_av_FRET_smooth)))*frame_rate/60);

        track_metrics{i}.mean_FRET = mean(trackfile{i}.cell_av_FRET_smooth,'omitnan');
        track_metrics{i}.std_FRET = std(trackfile{i}.cell_av_FRET_smooth,'omitnan');

        track_metrics{i}.avg_height = mean(track_metrics{i}.peak_metrics.peak_height);
        track_metrics{i}.avg_width = mean(track_metrics{i}.peak_metrics.peak_width);
        track_metrics{i}.avg_prominence = mean(track_metrics{i}.peak_metrics.peak_prominence);
        track_metrics{i}.avg_baseline = mean(track_metrics{i}.peak_metrics.peak_baseline);

        track_metrics_summary_table.track_nb(i) = i;
        track_metrics_summary_table.nb_of_peaks(i) = height(track_metrics{i}.peak_metrics);
        track_metrics_summary_table.peak_frequency(i) = track_metrics{i}.peak_frequency;
        track_metrics_summary_table.mean_FRET(i) = track_metrics{i}.mean_FRET;
        track_metrics_summary_table.std_FRET(i) = track_metrics{i}.std_FRET;
        track_metrics_summary_table.avg_height(i) = track_metrics{i}.avg_height;
        track_metrics_summary_table.avg_width(i) = track_metrics{i}.avg_width;
        track_metrics_summary_table.avg_prominence(i) = track_metrics{i}.avg_prominence;
        track_metrics_summary_table.avg_baseline(i) = track_metrics{i}.avg_baseline;
        track_metrics_summary_table.track_length(i) = max(size(trackfile{i}.cell_av_FRET_smooth));

    end
end

