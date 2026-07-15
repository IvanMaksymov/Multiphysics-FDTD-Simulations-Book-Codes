# vitis_hls -f build_hls.tcl

open_project fdtd_tmz_hls

set_top fdtd_tmz_step

add_files fdtd_tmz_hls.cpp
add_files -tb fdtd_tmz_tb.cpp


open_solution solution1

set_part xc7a100tfgg484-1

create_clock -period 10


csim_design

csynth_design

export_design -rtl verilog -format ip_catalog

exit


