# EPI_Stability_Phantom_QA
 Matlab code for performing QA in EPI data for EPI Stability
 
 Download all the Matlab files and execute the code "Main_Code_QA_EPIStability.m".
 Follow the instructions, which requires three inputs (at least) from the user:
 1. Select (via the Matlab uigetfolder() function) the PARENT DIRECTORY (i.e., the 
 folder that contains the three required measurement folders).
 2. Is the data being analysed Dicom or Nii?
 3. Do you want to use the defaults (i.e. everything runs automatically) or not?
 
 If the defaults are executed, four images will be created at the end:
 1. Slice views for each measurement
 2. A summary Figure for the four metrics: % fluctuation, drift, SFNR and Rdc
 
 If the defaults are not selected, then the code will require from you for each 
 measurement:
 1. If check ROI will be done or not
 2. If analysis results will be saved in the input or current work directory folder
 3. Repetion time if it is not accesible neither from the DICOM or Nii dataset
 4. Total number of slices if it is not accesible neither from the DICOM or Nii dataset
 5. Slices to be analysed if it is not accesible neither from the DICOM or Nii dataset
 6. Number of measurements to be analysed if it is not accesible neither from the DICOM or Nii dataset
 
 
