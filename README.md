# SSDS
Subject-specific diffusion segmentation pipeline 

Subject Specific Diffusion Statistics (SSDS)
adriana.azor16@imperial.ac.uk
Developed by Adriana Azor 
at C3NL, Imperial College London, London, UK
Last update: 20/03/2021
 
Description
This notebook wraps FSL & MRTRIX3 native commands in a pipeline to generate FA summary measures from whole tracts and boundary in native diffusion space

Dependencies
•	FSL installation
•	MRTRIX3 installation

Input data Raw DWI & T1w
Data organization:
  •	Working directory (with subjects.txt file list of all subjects in analysis)
      	Subject directory
          o	DWI (required)
          o	Raw diffusion file & bvec, bval (required)
          o	T1 (required)
          o	gre_map_field.nii.gz (if available)
          
          
          
 If fieldmap is available --> proceed with pipeline_withFM.sh
 If fieldmap is not available --> proceed with pipeline_noFM.sh
