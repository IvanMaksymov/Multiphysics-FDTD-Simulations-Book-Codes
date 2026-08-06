`timescale 1ns / 1ps

module fdtd_tmz
(
    input  wire       clk,
    input  wire       rst,
    output wire [5:0] LED
);

    parameter NX = 16;
    parameter NY = 16;
    parameter SIZE = NX*NY;

    //--------------------------------------------------
    // Fields
    //--------------------------------------------------

    reg signed [15:0] Ez [0:SIZE-1];
    reg signed [15:0] Hx [0:SIZE-1];
    reg signed [15:0] Hy [0:SIZE-1];

    integer i;

    //--------------------------------------------------
    // Coordinates
    //--------------------------------------------------

    reg [3:0] x;
    reg [3:0] y;

    wire [7:0] addr;
    wire [7:0] addr_r;
    wire [7:0] addr_u;

    assign addr   = y*NX + x;
    assign addr_r = y*NX + x + 1;
    assign addr_u = (y+1)*NX + x;

    //--------------------------------------------------
    // Fixed point coefficient
    // C = 0.25 in Q12
    //--------------------------------------------------

    localparam signed [15:0] C = 16'sd1024;

    reg signed [31:0] curl;

    //--------------------------------------------------
    // Clock divider
    // slows visual evolution
    //--------------------------------------------------

    reg [15:0] counter;

    wire step;

    assign step = (counter == 16'd5000);

    always @(posedge clk)
    begin
        if(rst)
            counter <= 0;

        else if(step)
            counter <= 0;

        else
            counter <= counter + 1;

    end

    //--------------------------------------------------
    // FDTD update
    //--------------------------------------------------

    always @(posedge clk)
    begin

        if(rst)
        begin

            x <= 1;
            y <= 1;

            for(i=0;i<SIZE;i=i+1)
            begin
                Ez[i] <= 0;
                Hx[i] <= 0;
                Hy[i] <= 0;
            end

        end

        else if(step)
        begin

            //--------------------------------------------------
            // update one cell each step
            //--------------------------------------------------
            if((x < NX-1) && (y < NY-1))
            begin

                // magnetic fields
                Hx[addr] <= Hx[addr]
                    - (((Ez[addr_u]-Ez[addr]) * C) >>> 12);


                Hy[addr] <= Hy[addr]
                    + (((Ez[addr_r]-Ez[addr]) * C) >>> 12);

                // electric field
                curl =
                    (Hy[addr]-Hy[addr-1])
                  - (Hx[addr]-Hx[addr-NX]);

                Ez[addr] <= Ez[addr]
                    + ((curl*C) >>> 12);

                //--------------------------------------------------
                // point source
                //--------------------------------------------------
                if((x==NX/2)&&(y==NY/2))
                    Ez[addr] <= Ez[addr] + 16'sd2000;

            end

            //--------------------------------------------------
            // move through grid
            //--------------------------------------------------
            if(x==NX-2)
            begin
                x <= 1;

                if(y==NY-2)
                    y <= 1;
                else
                    y <= y+1;

            end

            else
                x <= x+1;

        end

    end

    //--------------------------------------------------
    // LED display
    //--------------------------------------------------
    wire signed [15:0] centre;

    assign centre = Ez[(NY/2)*NX+(NX/2)];

    wire signed [15:0] abs_centre;

    assign abs_centre =
        centre[15] ? -centre : centre;

    assign LED[5:0] = abs_centre[15:10];


endmodule