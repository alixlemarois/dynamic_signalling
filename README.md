# dynamic_signalling
MATLAB code for ERK &amp; AKT biosensor quantification and dynamic pulse analysis.

The code in this repository enables to perform 3 steps of image analysis used in the Le Marois et al manuscript (https://www.biorxiv.org/content/10.1101/2024.05.14.594112v1). To run this code, an installation of Matlab 2022 or later is needed; however this can also be run using the Matlab online capability, accessible free of charge for those without a Matlab license – see here: https://uk.mathworks.com/products/matlab-online.html. Note that use time is limited to 20h/month.
The CODE folder contains 3 main modules and a demo script. Please refer to the header paragraph of each module for more detailed description.
1.	‘compute_biosensors’ extracts biosensor and morphological feature values from raw timelapse images with matched segmented masks.
2.	‘make_biosensor_tracks’ generates biosensor and morphological feature cell tracks from (x,y,t) track positions obtained from Trackmate: 
3.	Detection of ERK activity pulses and quantification of pulse features and pulse frequency: ‘detect_pulses’.
4.	The demo script can be run step by step to call the above modules and generate figures similar to those presented in Fig.3A of the manuscript: ‘demo_code’. The figures and data generated by this script are saved in the ‘DATA’ folder as each section is run.

Dependencies:

The bioformats-matlab integration is required to run this code and available here : https://www.openmicroscopy.org/bio-formats/. Add the 'bfmatlab' folder inside the 'CODE' folder for this to work properly.

The ‘notBoxPlot’ and ‘sigstar’ functions are required and included heres. For more details about these dependencies please refer to:
https://uk.mathworks.com/matlabcentral/fileexchange/39696-raacampbell-sigstar
https://uk.mathworks.com/matlabcentral/fileexchange/26508-notboxplot

Exemplar data is provided in the ‘DATA’ folder, and consists of 3 time-lapse images of H1975 cells expressing the EKAREV-NLS ERK-FRET sensor, with matched segmentation masks generated with ilastik. Each image has 2 channels – CFP and YFP, and 236 time points corresponding to approx. 15h of imaging. 
IMPORTANT NOTE: to enable easy running of this code, here no track correction is performed, as this would normally be performed manually by the user, to remove mitotic intervals and aberrant tracks. These steps were removed from this package for simplicity, meaning the pulse frequency quantification is different from the one performed in the manuscript, which contains additional manual correction steps. However, the quantification of the uncorrected data generated here does not differ significantly from what was obtained after track corrections.
