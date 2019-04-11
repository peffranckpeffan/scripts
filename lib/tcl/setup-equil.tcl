mol load psf ../common/solvate.psf pdb ../common/system.pdb
set all [atomselect top all]
$all set occupancy 0
$all set beta 0
$all writepdb refumb0.pdb

set all [atomselect top all]
$all set occupancy 0
$all set beta 0
set pr [atomselect top "segname A and backbone"]
$pr set occupancy 1.0
#measure center $pr weight mass
$all writepdb atoms.pdb

exit
