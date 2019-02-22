mol new ../../common/solvate.psf 
mol addfile ../../equil/out_eq2.restart.coor molid 0 waitfor all
mol addfile ../../SMD/out_smd.dcd molid 0 waitfor all
set range 20
set spacing 0.5
set totframes 5000

#intervalo em numeros de frames
set intv [expr double($totframes)/double($range)]
set med_intv [expr $intv/2]
set i 7
set b "75"
set fra [expr round([expr $i*$spacing*$intv + $spacing*$med_intv])]
set a [atomselect top all frame $fra]
$a writenamdbin out_smd-[format "%02d" [expr $i + 1]].1.restart.coor


#set win [expr int([expr $range/$spacing]) + 1] 
#numero de janelas
# set win [expr int([expr $range/$spacing]) + 1]    
#janela 00
# set a [atomselect top "segname B" frame 0]
# set center [measure center $a weight mass]
# $a writenamdbin out_smd-00.restart.coor

#  for {set i 0} {$i < $win} {incr i} {
#  	set j [format "%02d" [expr $i + 1]]
# # 	#comeca escrever a partir da janela 01
# 	set fra [expr round([expr $i*$spacing*$intv])]
# 	puts "$j $i $fra" 
# # 	set a [atomselect top all frame $fra]
# # 	$a writenamdbin out_smd-[format "%02d" $i + 1].restart.coor
# # 	puts $fra
#  } 

exit
