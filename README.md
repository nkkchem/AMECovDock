**AMECovDock - Automated Modeling Engine for covelent docking using ADFR**

This package holds everything needed to automate ADFR

In a covalently bound ligand, there are three atoms shared between the receptor and the ligand. It includes a common anchor atom in the receptor that connects with the ligand’s atoms. You can see both connected atoms in the PDBQT files of both the receptor and the ligand as well.

For with a receptor using  ADFR [1-2]  follows https://ccsb.scripps.edu/adfr/tutorial-covalent/

===============
**CURRENTLY SUPPORTED WARHEAD REACTIONS**
===============
acrylamides
alkynes
chloroacetamides
epoxyKetone_alphaSub
epoxyKetone_betaSub
esters
haloKetone_SN2
vinylmethyl ethers
vinylsulfones

==============
**OUTPUT**
==============
./Res_Dock/
    for ease of comparison and visualization of results
    receptorName_ligandName_reactionName_counter_jobName_out.pdbqt - the docked ligand
    scores.txt                                                     - list of docked ligands sorted by score

./output/
    additional information about docked ligands
    ./output/ligandName_reactionName_counter/
        .dro file          - from ADFR
        summary.dlg file   - from ADFR
        atoms.txt          - atom numbers chosen for docking, for additional validation

./pdb/
    for ease of re-docking and input validation
    .pdb    - result of AMECovDock_Reactions.py
    .pdbqt  - result of adfr prepare_ligand

===============
**ENVIRONMENT** - Constance
===============
module load python/anaconda3.2019.3
source /share/apps/python/anaconda3.2019.3/etc/profile.d/conda.sh
conda activate my-rdkit-env
module load gcc/6.1.0
module load openbabel/2.4.1

===============
DEPENDENCIES
===============
ADFR - http://adfr.scripps.edu/AutoDockFR/downloads.html
anaconda
rdkit - https://www.rdkit.org/docs/Install.html
gcc - used on Constance for compatibility with slurm

==============
NOTES - SETUP
==============
only 2 files need to be provided
  {ligand_file}.smi
    2 columns, no header: SMILES ID
  {receptor}_cov.trg
    for example run this command on your receptor.pdbqt (6wqf in this case)
    agfr -r 6wqf.pdbqt -b user -17.00 -5.00 15.00 40.00 40.00 40.00 -c 1379 1382 -t 1377 -x A:CYS145 -o 6wqf_cov

in AMECovDock_ADFR.sh
  modify  lines 14- 16     - receptor code, {ligand_file} base name, /path/to/ADFR/Installation/bin/
  modify lines 42-45       - Residue Number
in AMECovDock_Reactions.py
  modify line 99           - {ligand_file}.smi goes in supplier
  modify line 223          - SetResidueNumber

resulting files show
    receptorName - for keeping track of what receptorthe ligand was docked to
    ligandName   - for referencing ligand with original input
    reactionName - for easily identifying which warhead reacted with the receptor, may also be used to filter results according to interest
    counter      - when a reaction produces multiple unique products, the counter prevents them from overwriting

=================
TO RUN
./AMECovDock_ADFR.sh
=================

=================
References

1. Zhao, Y., Stoffler, D., & Sanner, M. (2006). Hierarchical and multi-resolution representation of protein flexibility. Bioinformatics, 22(22), 2768-2774.
2. Ravindranath, P. A., Forli, S., Goodsell, D. S., Olson, A. J., & Sanner, M. F. (2015). AutoDockFR: advances in protein-ligand docking with explicitly specified binding site flexibility. PLoS computational biology, 11(12), e1004586.
