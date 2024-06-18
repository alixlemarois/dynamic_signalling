
% Run this script section by section (use the 'run section' command in toolbar above while cursor is in one code section) to process raw data and obtain false
% colour FRET images, generate FRET tracks and detect ERK pulses.

%% SECTION 0 - initialise
addpath('Code')
addpath('Code/bfmatlab')
addpath('Data')

root_dir = pwd

%% SECTION I - Process raw images to obtain biosensor and morphological data + display example FRET images.
% this step may take a few minutes as the code runs through each image and
% each time point.
chdir('Data')
data_dir = pwd
compute_biosensors(pwd, 0, 0, 20, [50 50], 'ilastik', 2, 'C_T')

% Display FRET images
H1975DN_DMSO_FRET = load(strcat(data_dir, '/H1975DN_DMSO_segmented_props'), 'cell_av_FRET_image');
H1975DN_osi_FRET = load(strcat(data_dir, '/H1975DN_osi_segmented_props'), 'cell_av_FRET_image');
H1975DTP_osi_FRET = load(strcat(data_dir, '/H1975DTP_osi_segmented_props'), 'cell_av_FRET_image');

% t = 6h of imaging (90th time frame)
FRET_images = figure
subplot(1,3,1)
imshow(H1975DN_DMSO_FRET.cell_av_FRET_image(:,:,90))
caxis([0.6 1.6])
colormap(vertcat([0 0 0], jet(256)))
title('H1975-DN - DMSO')

subplot(1,3,2)
imshow(H1975DN_osi_FRET.cell_av_FRET_image(:,:,90))
caxis([0.6 1.6])
colormap(vertcat([0 0 0], jet(256)))
title('H1975-DN - osi')

subplot(1,3,3)
imshow(H1975DTP_osi_FRET.cell_av_FRET_image(:,:,90))
caxis([0.6 1.6])
colormap(vertcat([0 0 0], jet(256)))
title('H1975-DTP - osi')

% save image
saveas(FRET_images, 'FRET_images_6h.png')

%%  SECTION II - Import trackmate tracks and generate biosensor tracks + display example tracks

H1975DN_DMSO_trackmate_output = readmatrix('H1975DN_DMSO_trackmate.csv');
H1975DN_DMSO_trackmate_output = H1975DN_DMSO_trackmate_output(2:end, 3:9);

H1975DN_osi_trackmate_output = readmatrix('H1975DN_osi_trackmate.csv');
H1975DN_osi_trackmate_output = H1975DN_osi_trackmate_output(2:end, 3:9);

H1975DTP_osi_trackmate_output = readmatrix('H1975DTP_osi_trackmate.csv');
H1975DTP_osi_trackmate_output = H1975DTP_osi_trackmate_output(2:end, 3:9);

% import segmented features from previous step
H1975DN_DMSO = load(strcat(data_dir, '/H1975DN_DMSO_segmented_props'));
H1975DN_osi = load(strcat(data_dir, '/H1975DN_osi_segmented_props'));
H1975DTP_osi = load(strcat(data_dir, '/H1975DTP_osi_segmented_props'));

H1975DN_DMSO_tracks = make_biosensor_tracks(H1975DN_DMSO.im_segmented,H1975DN_DMSO.nuclei_props,H1975DN_DMSO_trackmate_output);
H1975DN_osi_tracks = make_biosensor_tracks(H1975DN_osi.im_segmented,H1975DN_osi.nuclei_props,H1975DN_osi_trackmate_output);
H1975DTP_osi_tracks = make_biosensor_tracks(H1975DTP_osi.im_segmented,H1975DTP_osi.nuclei_props,H1975DTP_osi_trackmate_output);

% display example tracks - track IDs chosen for their length

example_traces = figure
subplot(3,3,1)
plot(H1975DN_DMSO_tracks{2}.frame*4/60,H1975DN_DMSO_tracks{2}.cell_av_FRET)
ylim([0.7 1.6])
title('H1975-DN - DMSO','FontSize', 14,'FontWeight','normal')
subplot(3,3,4)
plot(H1975DN_DMSO_tracks{14}.frame*4/60,H1975DN_DMSO_tracks{14}.cell_av_FRET)
ylim([0.7 1.6])
ylabel('ERK activity','FontSize', 14)
subplot(3,3,7)
plot(H1975DN_DMSO_tracks{38}.frame*4/60,H1975DN_DMSO_tracks{38}.cell_av_FRET)
ylim([0.7 1.6])

subplot(3,3,2)
plot(H1975DN_osi_tracks{2}.frame*4/60,H1975DN_osi_tracks{2}.cell_av_FRET)
ylim([0.7 1.6])
title('H1975-DN - osi','FontSize', 14,'FontWeight','normal')
subplot(3,3,5)
plot(H1975DN_osi_tracks{17}.frame*4/60,H1975DN_osi_tracks{17}.cell_av_FRET)
ylim([0.7 1.6])
subplot(3,3,8)
plot(H1975DN_osi_tracks{34}.frame*4/60,H1975DN_osi_tracks{34}.cell_av_FRET)
ylim([0.7 1.6])
xlabel('time (h)','FontSize', 14)

subplot(3,3,3)
plot(H1975DTP_osi_tracks{4}.frame*4/60,H1975DTP_osi_tracks{4}.cell_av_FRET)
ylim([0.7 1.6])
title('H1975-DTP - osi','FontSize', 14,'FontWeight','normal')
subplot(3,3,6)
plot(H1975DTP_osi_tracks{16}.frame*4/60,H1975DTP_osi_tracks{16}.cell_av_FRET)
ylim([0.7 1.6])
subplot(3,3,9)
plot(H1975DTP_osi_tracks{34}.frame*4/60,H1975DTP_osi_tracks{34}.cell_av_FRET)
ylim([0.7 1.6])

%save data and images
saveas(example_traces, 'example_FRET_traces.png');
save('H1975_tracks.mat', 'H1975DN_DMSO_tracks', 'H1975DN_osi_tracks', 'H1975DTP_osi_tracks');

%% SECTION III - Detect pulses in tracks + display pulse frequency boxplots

% pulse detection code - 4 minute frame rate, 4-point smoothing ( = 16
% minute moving average), min peak prominence of 0.05
[H1975DN_DMSO_track_metrics, H1975DN_DMSO_track_metrics_summary_table] = detect_pulses(H1975DN_DMSO_tracks, 0, 4, 'trackmate', 4, 0.05)
[H1975DN_osi_track_metrics, H1975DN_osi_track_metrics_summary_table] = detect_pulses(H1975DN_osi_tracks, 0, 4, 'trackmate', 4, 0.05)
[H1975DTP_osi_track_metrics, H1975DTP_osi_track_metrics_summary_table] = detect_pulses(H1975DTP_osi_tracks, 0, 4, 'trackmate', 4, 0.05)

% display boxplots of pulse frequency.
fq_boxplots = figure
notBoxPlot(H1975DN_DMSO_track_metrics_summary_table.peak_frequency,0.5)
hold on
notBoxPlot(H1975DN_osi_track_metrics_summary_table.peak_frequency,1)
notBoxPlot(H1975DTP_osi_track_metrics_summary_table.peak_frequency,1.5)
ylim([-0.1 2.5])
xlim([0.25 1.75])
set(gca,'FontSize', 16)
ylabel('pulse frequency (h^{-1})')
set(gca,'XTick', [0.5 1 1.5],'XTickLabel', {'H1975-DN - DMSO','H1975-DN - osi','H1975-DTP - DMSO'})

% statistical tests
[~,p_DN_DMSO_v_DN_osi] = ttest2(H1975DN_DMSO_track_metrics_summary_table.peak_frequency, H1975DN_osi_track_metrics_summary_table.peak_frequency);
[~,p_DN_DMSO_v_DTP_osi] = ttest2(H1975DN_DMSO_track_metrics_summary_table.peak_frequency, H1975DTP_osi_track_metrics_summary_table.peak_frequency);
[~,p_DN_osi_v_DTP_osi] = ttest2(H1975DN_osi_track_metrics_summary_table.peak_frequency, H1975DTP_osi_track_metrics_summary_table.peak_frequency);

sigstar([0.5 0.95], p_DN_DMSO_v_DN_osi, [], 2.3)
sigstar([0.5 1.5], p_DN_DMSO_v_DTP_osi, [], 2.5)
sigstar([1.05 1.5], p_DN_osi_v_DTP_osi, [], 2.3)

% save data and images
saveas(fq_boxplots, 'pulse_frequency_boxplots.png')
save('H1975_pulse_analysis.mat', 'H1975DN_DMSO_track_metrics_summary_table', 'H1975DN_osi_track_metrics_summary_table', 'H1975DTP_osi_track_metrics_summary_table')
