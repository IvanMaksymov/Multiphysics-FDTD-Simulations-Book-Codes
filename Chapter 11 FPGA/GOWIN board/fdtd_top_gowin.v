`timescale 1ns / 1ps

module fdtd_top
(
    input  wire Clock,
    output wire IO_voltage
);

wire [5:0] LED;

wire busy;
wire done;

wire signed [15:0] ez_center;


// reset pulse after power-up
reg rst = 1'b1;
reg [15:0] rst_counter = 16'd0;


always @(posedge Clock)
begin
    if (rst_counter < 16'hffff)
    begin
        rst_counter <= rst_counter + 1'b1;
        rst <= 1'b1;
    end
    else
    begin
        rst <= 1'b0;
    end
end



fdtd_tmz
#(
    .NX(24),
    .NY(24),
    .DW(16)
)
u_fdtd
(
    .clk(Clock),
    .rst(rst),

    // free-running operation
    .start(1'b1),

    // Q6.10 coefficients
    .ce(16'sd256),
    .chx(16'sd256),
    .chy(16'sd256),

    // excitation
    .source(16'sd128),

    .busy(busy),
    .done(done),

    .LED(LED),

    .ez_center(ez_center)
);


// use LED0 pin as output indicator
assign IO_voltage = LED[0];


endmodule