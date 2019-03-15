mol new "../../../common/solvate.psf" 
mol addfile $env(center_file) molid 0 waitfor all
#set num_frames [molinfo top get numframes]
# set all [atomselect top $env(selection) frame $num_frames]
set all [atomselect top $env(selection)]
set center_mass [measure center $all weight mass]

set fileExists [file exist "center.tmp"]
if {$fileExists} {
	file delete -force ["center.tmp"]
}

set file [open "center.tmp" w]
puts $file $center_mass
close $file

exit
