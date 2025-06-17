from mdsam import phasetech2dir

datafolder='./Test_files'
phasetechdata=PhaseTech2DIRData(datafolder,'20250312')
phasetechdata.save(path_to_save)
morephasedata=PhaseTech2DIRData(path_to_some_hdf5_file)
