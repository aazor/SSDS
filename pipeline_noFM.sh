# Step 1. Pre-processing of diffusion data
 #The following cells preprocess raw diffusion data and fit tensor

 #Required dependency: FSL , MRtrix

 subj=`cat subject.txt`
 echo "eddy correcting"
 eddy_correct "$subj" "$subj"_ec 0; # change to match b0 from scanner 
 echo "create nodif" 
 for i in `cat subjects.txt`
 do eddy_correct "$i"/DWI/data.nii.gz "$i"/DWI/data_ec.nii.gz 0;
 fslroi "$i"/DWI/data_ec "$i"/DWI/nodif.nii.gz 0 1;

 # change first value to match b0 from scanner

 echo "fdt rotate" 
 fdt_rotate_bvecs "$i"/DWI/*.bvec "$i"/DWI/rotated_bvec "$i"/DWI/data_ec.ecclog;
 echo "bet" 
 bet "$i"/DWI/nodif.nii.gz "$i"/DWI/nodif_brain -f 0.2 -g 0 -m; 

 # QA brain extract 

 echo "dtifit" 
 dtifit --data="$i"/DWI/data.nii.gz --out="$i"/DWI/dti --mask="$i"/DWI/nodif_brain_mask bvecs="$i"/DWI/rotated_bvec --bvals="$i"/DWI/"$i".bval -w;
 done

 #Step 2. Setup paths for analysis & setup directories & pre-process data
 #The following cells setup the current environment for analysis and the directories required
 # | EDIT HERE: Project Directory, Tract Directory, Standard Directory. Export working directory, tract directory, and standard atlas directory.

 for i in `cat subjects.txt`;
 do mkdir ${workdir}/"$i"/boundary_mask;
 mkdir ${workdir}/"$i"/boundary_mask/MASK/;
 mkdir ${workdir}/"$i"/boundary_mask/transform;
 mkdir ${workdir}/"$i"/boundary_mask/tmp_files
 cp ${workdir}/"$i"/T1/T1.nii.gz ${workdir}/"$i"/boundary_mask/MASK
 cp ${workdir}/"$i"/DWI/nodif_brain.nii.gz ${workdir}/"$i"/boundary_mask/
 cp ${workdir}/"$i"/DWI/nodif.nii.gz ${workdir}/"$i"/boundary_mask/
 done

 echo " Directories ready"

 echo " **** REALIGNING BO IMAGES ***** "

 for i in `cat subjects.txt`;
 do flirt -dof 6 -in ${workdir}/"$i"/boundary_mask/nodif_brain.nii.gz -ref   ${workdir}/"$i"/DWI/data.nii.gz -out ${workdir}/"$i"/boundary_mask/nodif_brain.nii.gz -cost mutualinfo;
 done


 echo " **** Skull stripping of T1 ****
 	    -----------------------------------"

 for i in `cat subjects.txt`;
 do bet ${workdir}/"$i"/boundary_mask/MASK/T1.nii.gz ${workdir}/"$i"/boundary_mask/MASK/T1_brain.nii.gz -f 0.2 -B;
 done 

 #Step 3-b. Registration of T1 to B0 when field map is not available
 #The following cells runs the registration of T1 to B0 when field map is not available using BBR. This requires an initial segmentation of the T1 and estimation of a boundary mask that will then be used with a BBR cost function and a standard affine registration.

 
  echo  " **** Segmentation in T1 and preparation of WM maps for BBR ****
   ---------------------------------------------------------------------"

 for i in `cat subjects.txt`;
 do 5ttgen -nocrop fsl ${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz ${workdir}/"$i"/SSDS/MASK/5ttseg.mif -premasked -force;
 mrtransform -interp nearest -template ${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz ${workdir}/"$i"/SSDS/MASK/5ttseg.mif ${workdir}/"$i"/SSDS/MASK/5ttseg2.mif -force;
 mrconvert -coord 3 2 ${workdir}/"$i"/SSDS/MASK/5ttseg2.mif ${workdir}/"$i"/SSDS/MASK/wmseg.mif -force;
 mrconvert ${workdir}/"$i"/SSDS/MASK/wmseg.mif ${workdir}/"$i"/SSDS/MASK/wmseg.nii.gz -force
 done 

  #Boundary-based registration of T1 to B0

 echo  " **** Calculating matrix to move T1 to B0 using BBR ****
 	        -------------------------------------------------------"

 for i in `cat subjects.txt`;
 do 
 fslmaths ${workdir}/"$i"/SSDS/MASK/wmseg.nii.gz -bin ${workdir}/"$i"/SSDS/MASK/wmseg.nii.gz #binarize WM mask
 flirt -dof 6 -in ${workdir}/"$i"/SSDS/nodif_brain.nii.gz -ref ${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz -omat ${workdir}/"$i"/SSDS/transform/nodif2T1.mat -interp nearestneighbour; #estimate initial transform matrix from nodif to T1
 done

 for i in `cat subjects.txt`;
 do flirt -dof 6 -in ${workdir}/"$i"/SSDS/nodif_brain.nii.gz -ref ${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz -out ${workdir}/"$i"/SSDS/MASK/nodif2T1_brain.nii.gz -wmseg ${workdir}/"$i"/SSDS/MASK/wmseg.nii.gz -init ${workdir}/"$i"/SSDS/transform/nodif2T1.mat -omat ${workdir}/"$i"/SSDS/transform/nodif2T1_bbr.mat  -cost bbr -schedule /usr/local/fsl/etc/flirtsch/bbr.sch; 

 estimate transform matrix using BBR from nodif to T1

 convert_xfm -omat ${workdir}/"$i"/SSDS/transform/T1_to_nodif_bbr.mat -inverse ${workdir}/"$i"/SSDS/transform/nodif2T1_bbr.mat 

 #reverse matrix to obtain transformation of T1 to nodif

 flirt -in ${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz -ref ${workdir}/"$i"/SSDS/nodif_brain.nii.gz -out ${workdir}/"$i"/SSDS/MASK/T1_B0_bbr.nii.gz -applyxfm -init ${workdir}/"$i"/SSDS/transform/T1_to_nodif_bbr.mat 
 
 #apply BBR matrix to T1 image to move ot to nodif space and QC

 done


 #Step 4. Now working in diffusion space to estimate a 1-erosion whole-brain boundary map
 #The following cells runs the segmentation of the T1 brain image in diffusion space before estimating a whole brain boundary mask, a mask of the WM/CSF and WM/GM boundaries. Whole brain boundary is SSDS_B0. CSF boundary is SSDS_B0_csf, and GM boundary is SSDS_B0_csf. The last step exports the results in two separate CSV files.

 for i in `cat subjects.txt`;
 do 5ttgen -nocrop fsl ${workdir}/"$i"/SSDS/MASK/T1_B0_bbr.nii.gz ${workdir}/"$i"/SSDS/MASK/5ttseg_B0.mif -premasked -force;
 mrtransform -interp nearest -template ${workdir}/"$i"/SSDS/nodif_brain.nii.gz ${workdir}/"$i"/SSDS/MASK/5ttseg_B0.mif ${workdir}/"$i"/SSDS/MASK/5ttseg2_B0.mif -force;
 mrconvert -coord 3 2 ${workdir}/"$i"/SSDS/MASK/5ttseg2_B0.mif ${workdir}/"$i"/SSDS/MASK/wmseg_B0.mif -force;
 mrcalc ${workdir}/"$i"/SSDS/MASK/wmseg_B0.mif 0.5 -gt ${workdir}/"$i"/SSDS/MASK/wmseg_binary_B0.mif -force;
 maskfilter ${workdir}/"$i"/SSDS/MASK/wmseg_binary_B0.mif connect ${workdir}/"$i"/SSDS/MASK/wmseg_bin_connected_B0.mif -largest -force;
 maskfilter ${workdir}/"$i"/SSDS/MASK/wmseg_bin_connected_B0.mif erode ${workdir}/"$i"/SSDS/MASK/wmseg_bin_eroded_B0.mif -force;
 mrcalc ${workdir}/"$i"/SSDS/MASK/wmseg_bin_connected_B0.mif ${workdir}/"$i"/SSDS/MASK/wmseg_bin_eroded_B0.mif -sub ${workdir}/"$i"/SSDS/MASK/wmseg_bin_diff_B0.mif -force;
 mrconvert ${workdir}/"$i"/SSDS/MASK/wmseg_bin_diff_B0.mif ${workdir}/"$i"/SSDS/MASK/SSDS_B0.nii.gz -force;
 mrconvert ${workdir}/"$i"/SSDS/MASK/wmseg_B0.mif ${workdir}/"$i"/SSDS/MASK/wmseg_B0.nii.gz -force;
 fslmaths ${workdir}/"$i"/SSDS/MASK/wmseg_B0.nii.gz -thr 0.99 -bin ${workdir}/"$i"/SSDS/MASK/wmseg_B0.nii.gz ;
 done
 echo "BOUNDARY MASKS READY"

 echo "Estimating CSF and GM boundaries in WM map"

 for i in `cat subjects.txt`;
 do fslmaths ${workdir}/"$i"/SSDS/boundaries_2_B0/csfseg_B0.nii.gz -dilM ${workdir}/"$i"/SSDS/boundaries_2_B0/csfseg_B0_dil.nii.gz
 fslmaths ${workdir}/"$i"/SSDS/boundaries_2_B0/gmseg_B0.nii.gz -sub ${workdir}/"$i"/SSDS/boundaries_2_B0/csfseg_B0_dil.nii.gz ${workdir}/"$i"/SSDS/boundaries_2_B0/gmseg_B0_nocsf.nii.gz
 fslmaths ${workdir}/"$i"/SSDS/boundaries_2_B0/gmseg_B0_nocsf.nii.gz -thr 0 -uthr 1 -bin ${workdir}/"$i"/SSDS/boundaries_2_B0/gmseg_B0_nocsf.nii.gz
 fslmaths ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0.nii.gz -mas ${workdir}/"$i"/SSDS/boundaries_2_B0/csfseg_B0_dil.nii.gz ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_csf.nii.gz
 fslmaths ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0.nii.gz -sub ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_csf.nii.gz ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_gm.nii.gz
 fslmaths ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_gm.nii.gz -thr 0 -uthr 1 -bin ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_gm.nii.gz
 Done

 for i in `cat subjects.txt`;do stats=`fslstats ${i}/DWI/dti_FA.nii.gz -k ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_gm.nii.gz -M`;echo ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_gm.nii.gz ${i} ${stats};done>FA_results_GM_boundary.csv;
 for i in `cat subjects.txt`;do stats=`fslstats ${i}/DWI/dti_FA.nii.gz -k ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_csf.nii.gz -M`;echo ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_csf.nii.gz ${i} ${stats};done>FA_results_CSF_boundary.csv

 for i in `cat subjects.txt`;do stats=`fslstats ${i}/DWI/dti_FA.nii.gz -k ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0.nii.gz -M`;echo ${workdir}/"$i"/SSDS/boundaries_2_B0/SSDS_B0_csf.nii.gz ${i} ${stats};done>FA_results_boundary.csv

 #Step 5. Non-linear registration of MNI152 template to diffusion space
 #This step creates a non-linear warp to estimate the registration of the template to the native diffusion space. This template will later be used to move the tracts to the diffusion space.

 	echo "**** Create transformation matrix to move T1 to JHU - nonlinear ****
 	    -----------------------------------------------------------------------"

 for i in `cat subjects.txt`;
 do flirt -in ${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz -ref ${standard}/MNI152_T1_1mm_brain.nii.gz  -omat ${workdir}/"$i"/SSDS/transform/T1_2_MNI_aff_brain.mat -interp nearestneighbour;
 fnirt --in=${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz --ref=${standard}/MNI152_T1_1mm_brain.nii.gz  --iout=${workdir}/"$i"/SSDS/MASK/T1_2_MNI_aff_diff_brain.nii.gz --aff=${workdir}/"$i"/SSDS/transform/T1_2_MNI_aff_brain.mat --refmask=/Users/aa6616/Desktop/All_data/My_project/TEDS/Boundaries/MNI152_T1_1mm_brain_mask.nii.gz --cout=${workdir}/"$i"/SSDS/transform/T1_2_MNI_brain.nii.gz
 done

 echo "**** Inverse warp (from T1 --> MNI to MNI --> T1) ****
 	    -----------------------------------------------------------------------"

 for i in `cat subjects.txt`;
 do invwarp --ref=${workdir}/"$i"/SSDS/MASK/T1_brain.nii.gz --warp=${workdir}/"$i"/SSDS/transform/T1_2_MNI_brain.nii.gz --out=${workdir}/"$i"/SSDS/transform/MNI_2_T1_brain.nii.gz ;
 done


 # Step 6. Non-linear registration of tracts to diffusion space
 #This step uses the non-linear warp previously estimated to move the tracts from standard space to native diffusion space. The tracts used are initially eroded. Once registered, they are further thresholded and cross masked with the WM map to ensure exclusion of partial volume

  echo "Moving JHU tracts to native diffusion"

 # Create tract directory in main subject directories

 for i in `cat subjects.txt`;
 do 
 for f in `cat ${tractdir}/tracts.txt`; 
 do mkdir ${workdir}/"$i"/SSDS/tracts_ero/;
 cp ${tractdir}/"$f".nii.gz ${workdir}/"$i"/SSDS/tracts_ero/;
 done
 done

 #Warp tracts from MNI to B0 and clen tracts

 for i in `cat subjects.txt`;
 do  
 for f in `cat ${tractdir}/tracts.txt` 
 do applywarp --ref=${workdir}/"$i"/SSDS/nodif_brain.nii.gz --in=${tractdir}/"$f" --warp=${workdir}/"$i"/SSDS/transform/MNI_2_T1_brain.nii.gz --postmat=${workdir}/"$i"/SSDS/transform/T1_to_nodif_bbr.mat --out=${workdir}/"$i"/SSDS/boundaries_2_B0/tracts/"$f"_B0 --interp=nn;
 fslmaths ${workdir}/"$i"/SSDS/MASK/SSDS_B0.nii.gz -mul ${workdir}/"$i"/SSDS/tracts_ero/"$f"_B0.nii.gz ${workdir}/"$i"/SSDS/tracts_ero/"$f"_edge #cross-mask with boundary mask
 fslmaths ${workdir}/"$i"/SSDS/MASK/wmseg_B0.nii.gz -mul ${workdir}/"$i"/SSDS/tracts_ero/"$f"_B0 ${workdir}/"$i"/SSDS/tracts_ero/"$f"_B0_WM;
 fslmaths ${workdir}/"$i"/SSDS/tracts_ero/"$f"_B0_WM -thr 0.99 ${workdir}/"$i"/SSDS/tracts_ero/"$f"_B0_WM_0.99;
 fslmaths ${workdir}/"$i"/SSDS/tracts_ero/"$f"_B0_WM_0.99 -bin ${workdir}/"$i"/SSDS/tracts_ero/"$f"_center;
 done
 done 

 #Step 7. Cleaning up unwanted files

 echo "moving temporary files to tmp folder"
 for i in `cat subjects.txt`;
 do mv ${workdir}/"$i"/SSDS/tracts_ero/*T1_bin.nii.gz ${workdir}/"$i"/SSDS/tmp_files
 mv ${workdir}/"$i"/SSDS/tracts_ero/*warp.nii.gz ${workdir}/"$i"/SSDS/tmp_files
 mv ${workdir}/"$i"/SSDS/tracts_ero/*dil.nii.gz ${workdir}/"$i"/SSDS/tmp_files
 mv ${workdir}/"$i"/SSDS/MASK/*.mif ${workdir}/"$i"/SSDS/tmp_files
 mv ${workdir}/"$i"/SSDS/tracts_ero/*thr* ${workdir}/"$i"/SSDS/tmp_files
 mv ${workdir}/"$i"/SSDS/tracts_ero/*bin* ${workdir}/"$i"/SSDS/tmp_files
 mv ${workdir}/"$i"/SSDS/tracts_ero/*WM* ${workdir}/"$i"/SSDS/tmp_files
 mv ${workdir}/"$i"/SSDS/tracts_ero/*sm* ${workdir}/"$i"/SSDS/tmp_files
 done

 #Step 8. Pull data for statistics
 #This step pulls the mean FA from each estimated tract. 

 	echo " **********************************
          *     ALL DONE:                  *
 	       *    NOW EXTRACTING FA VALUES    *
 	       **********************************"


 for i in `cat subjects.txt`;
 do cd ${workdir}/"$i"/SSDS/tracts_ero/ && 
 ls *center.nii.gz>> tracts_center.txt
 cd ${workdir}
 done

 for i in `cat subjects.txt`;do for roi in `cat ${workdir}/"$i"/SSDS/tracts_ero/tracts_center.txt`;do stats=`fslstats ${i}/DWI/dti_FA.nii.gz -k ${workdir}/"$i"/SSDS/tracts_ero/${roi} -M -S`;echo ${workdir}/"$i"/SSDS/tracts_ero/${roi} ${i} ${stats};done;done>FA_results.csv;

 echo "   *     ~~~~~~~~~~~~~~~~~~~~~~~COMPLETED~~~~~~~~~~~~~~~~~~~~~~~        *            "
 done
