mol new ../common/solvate.psf 
mol addfile ../equil/out_eq2.restart.coor molid 0 waitfor all
mol addfile ../SMD/out_smd.dcd molid 0 waitfor all
set range 20
set spacing 0.5
set totframes 5000
set intv [expr double($totframes)/double($range)] 
set win [expr int([expr $range/$spacing]) + 1] 
set a [atomselect top all frame 0]
$a writenamdbin out_smd-00.restart.coor
for {set i 0} {$i < $win} {incr i} {
set j [format "%02d" [expr $i + 1]]
set fra [expr round([expr $i*$spacing*$intv])]
set a [atomselect top all frame $fra]
$a writenamdbin out_smd-$j.restart.coor
puts $fra
} 

exit
