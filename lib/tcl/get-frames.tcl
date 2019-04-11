mol new ../common/system.psf 
mol addfile ../equil/out_eq2.restart.coor molid 0 waitfor all
mol addfile ../SMD/out_smd.dcd molid 0 waitfor all

set window_path $env(window_path)
set numframes [molinfo top get numframes]
set range $env(totrange)

set positions [split $env(window_path) ,]


set position_frame [atomselect top all frame 0]
$position_frame writenamdbin out_smd-00.restart.coor

set i 1

foreach position $positions {
	set spacing [expr $position - double([lindex $positions 0])]
	set fra [expr round([expr ($numframes*$spacing)/$range])]
	set position_frame [atomselect top all frame $fra]
	$position_frame writenamdbin out_smd-[format "%02d" [expr $i]].restart.coor
	set i [expr $i + 1]
}

 exit
