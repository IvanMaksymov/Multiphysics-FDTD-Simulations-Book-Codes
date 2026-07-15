set PROJECT fdtd_system
set HLS_IP_PATH "/media/ivan/Untitled/Ivan_Work_Dec2025/for_FDTD_book/FPGA/fdtd_tmz_hls/solution1/impl/ip"

create_project $PROJECT ./$PROJECT -part xc7a100tifgg484-1L -force

set_property ip_repo_paths [list $HLS_IP_PATH] [current_project]
update_ip_catalog

create_bd_design fdtd_bd

create_bd_cell -type ip -vlnv xilinx.com:hls:fdtd_tmz_step:1.0 fdtd_tmz_step_0

# Clock + reset
create_bd_port -dir I -type clk clk
create_bd_port -dir I -type rst resetn

connect_bd_net [get_bd_ports clk]    [get_bd_pins fdtd_tmz_step_0/ap_clk]
connect_bd_net [get_bd_ports resetn] [get_bd_pins fdtd_tmz_step_0/ap_rst_n]

# AXI-Lite external
make_bd_intf_pins_external [get_bd_intf_pins fdtd_tmz_step_0/s_axi_control]

validate_bd_design
save_bd_design

make_wrapper -files [get_files fdtd_bd.bd] -top
add_files -norecurse ./fdtd_system/fdtd_system.gen/sources_1/bd/fdtd_bd/hdl/fdtd_bd_wrapper.v
set_property top fdtd_bd_wrapper [current_fileset]

generate_target all [get_files fdtd_bd.bd]

launch_runs synth_1 -jobs 4
wait_on_run synth_1

launch_runs impl_1 -jobs 4
wait_on_run impl_1

open_run impl_1

# Set config voltage + relax IO DRCs
set_property CFGBVS VCCO [current_design]
set_property CONFIG_VOLTAGE 3.3 [current_design]
set_property SEVERITY {Warning} [get_drc_checks NSTD-1]
set_property SEVERITY {Warning} [get_drc_checks UCIO-1]

write_bitstream -force fdtd_bd_wrapper.bit

puts "=============================================="
puts " FDTD FPGA BITSTREAM GENERATED SUCCESSFULLY "
puts "=============================================="

