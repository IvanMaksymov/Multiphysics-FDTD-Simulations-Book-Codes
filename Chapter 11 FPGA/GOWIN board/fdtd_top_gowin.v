`timescale 1ns / 1ps

module top
(
    input wire clk,
    output wire [5:0] leds
);


fdtd_tmz u_fdtd
(
    .clk(clk),
    .rst(1'b0),
    .LED(leds)
);


endmodule