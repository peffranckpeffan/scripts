mol new ../common/solvate.psf 
mol addfile ../equil/out_eq2.restart.coor molid 0

set all [atomselect top all]
$all set beta 0.0
$all set occupancy 0.0
$all writepdb refumb0.pdb
set prot [atomselect top "segname A and backbone"]
$prot set occupancy 1.0
set protr [atomselect top "segname A and noh and same residue as within 6 of segname B"]
$protr set beta 3.0
measure center $prot weight mass
set am [atomselect top "segname B and backbone"]
$am set occupancy 2.0
measure center $am weight mass
set amr [atomselect top "segname B and noh"]
$amr set beta 4.0
$all writepdb atoms.pdb

exit
