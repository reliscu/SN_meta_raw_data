### README for raw FASTQ data used in Tran-Maynard, et al. Neuron 2021
### https://doi.org/10.1016/j.neuron.2021.09.001

- Data for the 14 samples as analyzed and presented in the preprint for this project (https://doi.org/10.1101/2020.10.07.329839) are here,
in their sample-specific directories, following "donorID_brainRegion[_neun]" [if that sample went through NeuN enrichment via FACS].
   # Sample Br5182_NAc has a '_reseq' for "re-sequenced", as it was found these was unintentionally adaptor-trimmed in the prep. of the FASTQ's
   #    (this was thus re-sequenced & corrected in the revision process)

- 10 additional samples were added in the revision, which are similarly organized in the Feb2021/ subdirectory.









## MNT - Internal notes for re-organization =====
# Move the index files (supplement from `bcl2fastq` output, which we didn't always get from the sequencing core) to a misc_or_QC/
#mkdir misc_or_QC/
#mv *_I* misc_or_QC/

# Remove any redundant copies of FASTQ files, already organized into sample subdirs (see AEJ's '../move_reads.R' 11Oct2019)
#rm 5161_AMY*
#rm 5161_DLPFC*
#rm 5161_HPC_suc*
#rm 5161_sACC*
#rm 5207_NAC*
#rm 5212_Amy*
#rm 5212_HPC*
#rm 5212_sACC*
#rm 5287_DLPFC_suc*
#rm 5287_NAC*
#rm BR5161DPC*	# (this was a typo by the core--see folder Br5161_HPC/)
#rm BR5287DLPFC*
#rm BR5161_Nac*
#rm BR5212_DLPFC*
#rm BR5212_Nac*
#rm BR5287_HPC*

# Move sucrose-gradient-processed samples, QC-removed, and Br5276_Amy (which was done with Br1469_Hb for QC'ing the in-house S3e FACS)
#							   (this sample was re-processed in Feb2021, and with NeuN enrichment)
#mv Br5276_Amy/ misc_or_QC/
#mv *_suc misc_or_QC/
#mv Br5287_DLPFC/ misc_or_QC/

# And finally Br5182_NAc/ (see above note for explanation)
#mv Br5182_NAc/ misc_or_QC/


