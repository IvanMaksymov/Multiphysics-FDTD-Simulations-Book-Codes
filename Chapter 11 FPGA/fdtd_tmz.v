`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// 2‑D TMz FDTD Solver (Educational RTL)
// - One cell updated per clock cycle
// - Fixed‑point arithmetic (signed Q6.10)
//////////////////////////////////////////////////////////////////////////////////

module fdtd_tmz
#(
    parameter NX = 24,
    parameter NY = 24,
    parameter DW = 16
)
(
    input  wire                 clk,
    input  wire                 rst,

    // Start one FDTD timestep
    input  wire                 start,

    // FDTD coefficients (Q6.10)
    input  wire signed [DW-1:0] ce,
    input  wire signed [DW-1:0] chx,
    input  wire signed [DW-1:0] chy,

    // Soft source excitation
    input  wire signed [DW-1:0] source,

    // Status
    output reg                  busy,
    output reg                  done,

    // Ez centre value
    output wire signed [DW-1:0] ez_center,

    // History read port
    input  wire [9:0]           hist_addr,
    output wire signed [DW-1:0] hist_data
);

    localparam SIZE = NX * NY;

    //////////////////////////////////////////////////////////////
    // Field memories (block RAM)
    //////////////////////////////////////////////////////////////
    (* ram_style = "block" *) reg signed [DW-1:0] Ez [0:SIZE-1];
    (* ram_style = "block" *) reg signed [DW-1:0] Hx [0:SIZE-1];
    (* ram_style = "block" *) reg signed [DW-1:0] Hy [0:SIZE-1];

    assign ez_center = Ez[(NY/2)*NX + (NX/2)];

    //////////////////////////////////////////////////////////////
    // Ez history buffer (1024 samples)
    //////////////////////////////////////////////////////////////
    (* ram_style = "block" *) reg signed [DW-1:0] ez_history [0:1023];
    reg [9:0] ez_hist_wrptr = 10'd0;

    assign hist_data = ez_history[hist_addr];

    //////////////////////////////////////////////////////////////
    // Controller states
    //////////////////////////////////////////////////////////////
    localparam S_CLEAR    = 3'd0;
    localparam S_IDLE     = 3'd1;
    localparam S_UPDATEHX = 3'd2;
    localparam S_UPDATEHY = 3'd3;
    localparam S_UPDATEEZ = 3'd4;
    localparam S_BOUNDARY = 3'd5;
    localparam S_DONE     = 3'd6;

    reg [2:0] state;

    //////////////////////////////////////////////////////////////
    // Grid coordinates
    //////////////////////////////////////////////////////////////
    reg [7:0] x;
    reg [7:0] y;

    //////////////////////////////////////////////////////////////
    // Address calculation
    //////////////////////////////////////////////////////////////
    wire [15:0] addr       = y * NX + x;
    wire [15:0] addr_up    = (y + 1) * NX + x;
    wire [15:0] addr_right = y * NX + (x + 1);
    wire [15:0] addr_left  = y * NX + (x - 1);
    wire [15:0] addr_ym1   = (y - 1) * NX + x;

    //////////////////////////////////////////////////////////////
    // Arithmetic registers
    //////////////////////////////////////////////////////////////
    reg signed [31:0] delta;
    reg signed [31:0] curl;

    //////////////////////////////////////////////////////////////
    // Main controller
    //////////////////////////////////////////////////////////////
    integer i;

    always @(posedge clk) begin
        if (rst) begin
            state <= S_CLEAR;
            busy  <= 1'b1;
            done  <= 1'b0;
            x     <= 0;
            y     <= 0;
            ez_hist_wrptr <= 0;

            // Clear memories
            for (i = 0; i < SIZE; i = i + 1) begin
                Ez[i] <= 0;
                Hx[i] <= 0;
                Hy[i] <= 0;
            end
        end
        else begin
            case (state)

                //////////////////////////////////////////////////////
                // Clear field memories
                //////////////////////////////////////////////////////
                S_CLEAR: begin
                    Ez[addr] <= 0;
                    Hx[addr] <= 0;
                    Hy[addr] <= 0;

                    if (x == NX-1) begin
                        x <= 0;
                        if (y == NY-1) begin
                            y    <= 0;
                            busy <= 1'b0;
                            state <= S_IDLE;
                        end else begin
                            y <= y + 1;
                        end
                    end else begin
                        x <= x + 1;
                    end
                end

                //////////////////////////////////////////////////////
                // Wait for start
                //////////////////////////////////////////////////////
                S_IDLE: begin
                    done <= 1'b0;

                    if (start) begin
                        busy  <= 1'b1;
                        x     <= 0;
                        y     <= 0;
                        state <= S_UPDATEHX;
                    end
                end

                //////////////////////////////////////////////////////
                // Update Hx
                //////////////////////////////////////////////////////
                S_UPDATEHX: begin
                    delta = Ez[addr_up] - Ez[addr];
                    Hx[addr] <= Hx[addr] - ((delta * chx) >>> 10);

                    if (y == NY-2) begin
                        y <= 0;
                        if (x == NX-1) begin
                            x     <= 0;
                            state <= S_UPDATEHY;
                        end else begin
                            x <= x + 1;
                        end
                    end else begin
                        y <= y + 1;
                    end
                end

                //////////////////////////////////////////////////////
                // Update Hy
                //////////////////////////////////////////////////////
                S_UPDATEHY: begin
                    delta = Ez[addr_right] - Ez[addr];
                    Hy[addr] <= Hy[addr] + ((delta * chy) >>> 10);

                    if (x == NX-2) begin
                        x <= 0;
                        if (y == NY-1) begin
                            x     <= 1;
                            y     <= 1;
                            state <= S_UPDATEEZ;
                        end else begin
                            y <= y + 1;
                        end
                    end else begin
                        x <= x + 1;
                    end
                end

                //////////////////////////////////////////////////////
                // Update Ez
                //////////////////////////////////////////////////////
                S_UPDATEEZ: begin
                    curl =
                        (Hy[addr] - Hy[addr_left]) -
                        (Hx[addr] - Hx[addr_ym1]);

                    Ez[addr] <= Ez[addr] + ((curl * ce) >>> 10);

                    if ((x == NX/2) && (y == NY/2)) begin
                        Ez[addr] <= Ez[addr] + ((curl * ce) >>> 10) + source;
                    end

                    if (x == NX-1) begin
                        x <= 1;
                        if (y == NY-1) begin
                            y     <= 0;
                            state <= S_BOUNDARY;
                        end else begin
                            y <= y + 1;
                        end
                    end else begin
                        x <= x + 1;
                    end
                end

                //////////////////////////////////////////////////////
                // PEC boundaries
                //////////////////////////////////////////////////////
                S_BOUNDARY: begin
                    Ez[y*NX]           <= 0;
                    Ez[y*NX + NX-1]    <= 0;
                    Ez[x]              <= 0;
                    Ez[(NY-1)*NX + x]  <= 0;

                    if (y == NY-1) begin
                        y     <= 0;
                        state <= S_DONE;
                    end else begin
                        y <= y + 1;
                    end
                end

                //////////////////////////////////////////////////////
                // Done
                //////////////////////////////////////////////////////
                S_DONE: begin
                    busy <= 1'b0;
                    done <= 1'b1;

                    // Store Ez centre sample
                    ez_history[ez_hist_wrptr] <= ez_center;
                    ez_hist_wrptr <= ez_hist_wrptr + 1'b1;

                    state <= S_IDLE;
                end

            endcase
        end
    end

endmodule

