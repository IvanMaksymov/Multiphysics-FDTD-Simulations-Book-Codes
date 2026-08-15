`timescale 1ns / 1ps

//=========================================================================
// Ternary FDTD (TMz), hardware port of the C reference model:
//
//   aHx[i][j] += Ez[i][j]   - Ez[i][j+1];
//   aHy[i][j] += Ez[i+1][j] - Ez[i][j];
//   Hx/Hy[i][j] = quantise(&aHx/aHy[i][j], old);
//
//   aEz[i][j] += (Hy[i][j]-Hy[i-1][j]) - (Hx[i][j]-Hx[i][j-1]);
//   Ez[i][j]  = quantise(&aEz[i][j], old);
//
//   aEz[sx][sy] += source(n);
//   Ez[sx][sy]  = quantise(&aEz[sx][sy], old);   <- quantised AGAIN
//
// quantise(): residual-preserving, hysteretic (holds old value when
// the accumulator is inside the +/-THRESHOLD band), matching the C
// exactly - not a stateless round.
//
// There is no Courant-coefficient multiply anywhere in this design -
// the C model doesn't have one, so neither does this hardware. Fields
// are plain ternary (-1,0,+1); the physics lives entirely in the
// integer accumulate + threshold-crossing behaviour.
//=========================================================================

module fdtd_tmz
(
    input  wire               clk,
    input  wire               rst,

    input  wire [12:0]        display_addr,   // WIDER: SIZE=6400 needs 13 bits
    output wire signed [15:0] display_ez,      // sign-extended from ternary storage

    output wire [5:0]         LED
);

parameter NX   = 64;
parameter NY   = 64;
parameter SIZE = NX * NY;   // 64*64

//=====================================================================
// FDTD PARAMETERS
//=====================================================================

parameter STEP_DELAY = 26'd20000000;

// Fixed quantiser threshold - matches C's #define THRESHOLD 4
localparam signed [15:0] THRESHOLD = 16'sd4;

// Source cell (matches C's sx = NX/2, sy = NY/2)
localparam [6:0] SX = NX/2;   // 40
localparam [6:0] SY = NY/2;   // 40

// Source: amplitude 4, frequency ~0.01 cycles/step (matches C's
// SOURCE_AMPLITUDE=4.0f, SOURCE_FREQUENCY=0.01f). Implemented as a
// 16-bit phase accumulator + 16-entry sine LUT (top 4 bits index the
// table) instead of a true sinf() call. PHASE_INCREMENT=655 gives
// 655/65536 = 0.009995 cycles/step, within ~0.05% of the C target -
// tune this constant if your toolchain needs a different match.
localparam [15:0] PHASE_INCREMENT = 16'd655;

//=====================================================================
// FIELD MEMORIES
//=====================================================================

// Ternary field state (-1, 0, +1). 8 bits used for readability/ease of
// debug even though 2 bits would technically suffice.
(* ram_style = "block" *)
reg signed [7:0] Ez_mem [0:SIZE-1];
(* ram_style = "block" *)
reg signed [7:0] Hx_mem [0:SIZE-1];
(* ram_style = "block" *)
reg signed [7:0] Hy_mem [0:SIZE-1];

// Residual accumulators - MUST persist in BRAM between timesteps,
// unlike the old 16-bit fixed-point design which had no separate
// residual concept. This is what makes the hysteretic quantiser work.
(* ram_style = "block" *)
reg signed [15:0] aEz_mem [0:SIZE-1];
(* ram_style = "block" *)
reg signed [15:0] aHx_mem [0:SIZE-1];
(* ram_style = "block" *)
reg signed [15:0] aHy_mem [0:SIZE-1];

//=====================================================================
// DISPLAY READ
//=====================================================================

reg [12:0] display_addr_reg;
reg signed [7:0] display_ez_reg;

always @(posedge clk) begin
    display_addr_reg <= display_addr;
    display_ez_reg   <= Ez_mem[display_addr_reg];
end

// Sign-extend the ternary value (-1,0,+1) to the 16-bit output port
// so existing downstream consumers (VGA driver etc.) don't need to
// change their bus width, even though the dynamic range is now tiny.
assign display_ez = {{8{display_ez_reg[7]}}, display_ez_reg};

//=====================================================================
// CLOCK DIVIDER
//=====================================================================

reg [25:0] counter;
wire step = (counter == STEP_DELAY);

always @(posedge clk) begin
    if (rst)
        counter <= 26'd0;
    else if (step)
        counter <= 26'd0;
    else
        counter <= counter + 1'b1;
end

//=====================================================================
// SOURCE
//=====================================================================

reg [15:0] source_phase;
reg signed [7:0] sine_value;

always @(posedge clk) begin
    if (rst)
        source_phase <= 16'd0;
    else if (step)
        source_phase <= source_phase + PHASE_INCREMENT;
end

// trunc(4*sin(2*pi*k/16)) for k=0..15 - matches C's (int)(4.0f*sinf(...))
always @* begin
    case (source_phase[15:12])
        4'd0:  sine_value =  8'sd0;
        4'd1:  sine_value =  8'sd1;
        4'd2:  sine_value =  8'sd2;
        4'd3:  sine_value =  8'sd3;
        4'd4:  sine_value =  8'sd4;
        4'd5:  sine_value =  8'sd3;
        4'd6:  sine_value =  8'sd2;
        4'd7:  sine_value =  8'sd1;
        4'd8:  sine_value =  8'sd0;
        4'd9:  sine_value = -8'sd1;
        4'd10: sine_value = -8'sd2;
        4'd11: sine_value = -8'sd3;
        4'd12: sine_value = -8'sd4;
        4'd13: sine_value = -8'sd3;
        4'd14: sine_value = -8'sd2;
        4'd15: sine_value = -8'sd1;
        default: sine_value = 8'sd0;
    endcase
end

//=====================================================================
// FDTD CONTROLLER STATES
//=====================================================================

localparam [4:0]
    S_INIT_CLEAR = 5'd0,   // zero every memory at reset, like C's static init
    S_IDLE       = 5'd1,
    S_H_READ     = 5'd2,
    S_H_WAIT     = 5'd3,
    S_H_CALC     = 5'd4,
    S_H_WRITE    = 5'd5,
    S_E_READ_H   = 5'd6,
    S_E_WAIT_H   = 5'd7,
    S_E_READ_EZ  = 5'd8,
    S_E_WAIT_EZ  = 5'd9,
    S_E_CALC     = 5'd10,
    S_E_SOURCE   = 5'd11,
    S_E_WRITE    = 5'd12,
    S_NEXT       = 5'd13,
    S_WAIT       = 5'd14;

reg [4:0] state;
reg [6:0] x, y;           // 0..79, needs 7 bits
reg h_phase;

reg [12:0] init_idx;

//=====================================================================
// ADDRESSES
//
// H phase reads/writes cell (x,y) and needs neighbours (x+1,y) and
// (x,y+1) - safe without clamping because the H sweep is bounded to
// x,y in [0, NX-2], so x+1 and y+1 never exceed NX-1/NY-1.
//
// E phase reads/writes cell (x,y) and needs neighbours (x-1,y) and
// (x,y-1) - safe without clamping because the E sweep is bounded to
// x,y in [1, NX-2].
//=====================================================================

wire [12:0] addr_cur = y * NX + x;
wire [12:0] addr_r   = y * NX + (x + 7'd1);   // H phase only
wire [12:0] addr_u   = (y + 7'd1) * NX + x;   // H phase only
wire [12:0] addr_l   = y * NX + (x - 7'd1);   // E phase only
wire [12:0] addr_d   = (y - 7'd1) * NX + x;   // E phase only

//=====================================================================
// WORKING REGISTERS
//=====================================================================

reg signed [7:0]  Ez_cur, Ez_r, Ez_u;
reg signed [7:0]  Hx_cur, Hx_d;
reg signed [7:0]  Hy_cur, Hy_l;
reg signed [15:0] aHx_cur, aHy_cur, aEz_cur;

reg signed [15:0] diff_x, diff_y, diff_e;
reg signed [15:0] acc_x,  acc_y,  acc_e, acc_e2;

reg signed [7:0]  Hx_next, Hy_next;
reg signed [15:0] aHx_next, aHy_next;

reg signed [7:0]  Ez_partial, Ez_final;
reg signed [15:0] aEz_partial, aEz_final;

//=====================================================================
// WRITE CONTROL
//=====================================================================

reg [12:0] write_addr;
reg signed [7:0]  write_data_ez, write_data_hx, write_data_hy;
reg signed [15:0] acc_data_ez,  acc_data_hx,  acc_data_hy;
reg write_enable;
reg write_to_h;
reg clearing;   // write_enable used for both init-clear and normal writes

//=====================================================================
// BRAM WRITE - single process for field + its residual accumulator
//=====================================================================

always @(posedge clk) begin
    if (write_enable) begin
        if (clearing) begin
            Ez_mem[write_addr]  <= 8'sd0;
            Hx_mem[write_addr]  <= 8'sd0;
            Hy_mem[write_addr]  <= 8'sd0;
            aEz_mem[write_addr] <= 16'sd0;
            aHx_mem[write_addr] <= 16'sd0;
            aHy_mem[write_addr] <= 16'sd0;
        end else if (write_to_h) begin
            Hx_mem[write_addr]  <= write_data_hx;
            Hy_mem[write_addr]  <= write_data_hy;
            aHx_mem[write_addr] <= acc_data_hx;
            aHy_mem[write_addr] <= acc_data_hy;
        end else begin
            Ez_mem[write_addr]  <= write_data_ez;
            aEz_mem[write_addr] <= acc_data_ez;
        end
    end
end

//=====================================================================
// FDTD STATE MACHINE
//=====================================================================

always @(posedge clk) begin
    if (rst) begin
        state        <= S_INIT_CLEAR;
        init_idx     <= 13'd0;
        x            <= 7'd0;
        y            <= 7'd0;
        h_phase      <= 1'b1;
        write_enable <= 1'b0;
        clearing     <= 1'b1;
    end else begin
        write_enable <= 1'b0;

        case (state)

        //=============================================================
        // POWER-ON / RESET CLEAR - zero every memory, matching the
        // C model's zero-initialised static arrays. Also guarantees
        // the never-swept boundary cells read back as 0 forever,
        // which is how the C model's Ez[i][0]=0 etc. boundary is
        // actually achieved (those cells are simply never written by
        // the interior sweeps either way).
        //=============================================================
        S_INIT_CLEAR: begin
            clearing     <= 1'b1;
            write_addr   <= init_idx;
            write_enable <= 1'b1;

            if (init_idx == SIZE-1) begin
                init_idx <= 13'd0;
                state    <= S_IDLE;
            end else begin
                init_idx <= init_idx + 13'd1;
            end
        end

        S_IDLE: begin
            clearing <= 1'b0;
            if (step) begin
                x       <= 7'd0;
                y       <= 7'd0;
                h_phase <= 1'b1;
                state   <= S_H_READ;
            end
        end

        //=============================================================
        // H FIELD UPDATE - sweep x,y over [0, NX-2] (matches C's
        // "for(i=0;i<NX-1;i++)"), no Courant multiply, just raw
        // ternary neighbour differences.
        //=============================================================
        S_H_READ: begin
            state <= S_H_WAIT;
        end

        S_H_WAIT: begin
            Ez_cur  <= Ez_mem[addr_cur];
            Ez_r    <= Ez_mem[addr_r];
            Ez_u    <= Ez_mem[addr_u];
            Hx_cur  <= Hx_mem[addr_cur];
            Hy_cur  <= Hy_mem[addr_cur];
            aHx_cur <= aHx_mem[addr_cur];
            aHy_cur <= aHy_mem[addr_cur];
            state   <= S_H_CALC;
        end

        S_H_CALC: begin
            // aHx += Ez[i][j] - Ez[i][j+1]  ->  diff_x = Ez_cur - Ez_u
            // aHy += Ez[i+1][j] - Ez[i][j]  ->  diff_y = Ez_r - Ez_cur
            diff_x = Ez_cur - Ez_u;
            diff_y = Ez_r   - Ez_cur;
            acc_x  = aHx_cur + diff_x;
            acc_y  = aHy_cur + diff_y;

            if (acc_x >= THRESHOLD) begin
                Hx_next  = 8'sd1;
                aHx_next = acc_x - THRESHOLD;
            end else if (acc_x <= -THRESHOLD) begin
                Hx_next  = -8'sd1;
                aHx_next = acc_x + THRESHOLD;
            end else begin
                Hx_next  = Hx_cur;
                aHx_next = acc_x;
            end

            if (acc_y >= THRESHOLD) begin
                Hy_next  = 8'sd1;
                aHy_next = acc_y - THRESHOLD;
            end else if (acc_y <= -THRESHOLD) begin
                Hy_next  = -8'sd1;
                aHy_next = acc_y + THRESHOLD;
            end else begin
                Hy_next  = Hy_cur;
                aHy_next = acc_y;
            end

            state <= S_H_WRITE;
        end

        S_H_WRITE: begin
            write_addr    <= addr_cur;
            write_data_hx <= Hx_next;
            write_data_hy <= Hy_next;
            acc_data_hx   <= aHx_next;
            acc_data_hy   <= aHy_next;
            write_to_h    <= 1'b1;
            write_enable  <= 1'b1;
            state         <= S_NEXT;
        end

        //=============================================================
        // E FIELD UPDATE - sweep x,y over [1, NX-2] (matches C's
        // "for(i=1;i<NX-1;i++)").
        //=============================================================
        S_E_READ_H: begin
            state <= S_E_WAIT_H;
        end

        S_E_WAIT_H: begin
            Hy_cur <= Hy_mem[addr_cur];
            Hy_l   <= Hy_mem[addr_l];
            Hx_cur <= Hx_mem[addr_cur];
            Hx_d   <= Hx_mem[addr_d];
            state  <= S_E_READ_EZ;
        end

        S_E_READ_EZ: begin
            state <= S_E_WAIT_EZ;
        end

        S_E_WAIT_EZ: begin
            Ez_cur  <= Ez_mem[addr_cur];
            aEz_cur <= aEz_mem[addr_cur];
            state   <= S_E_CALC;
        end

        S_E_CALC: begin
            // aEz += (Hy[i][j]-Hy[i-1][j]) - (Hx[i][j]-Hx[i][j-1])
            diff_e = (Hy_cur - Hy_l) - (Hx_cur - Hx_d);
            acc_e  = aEz_cur + diff_e;

            if (acc_e >= THRESHOLD) begin
                Ez_partial  = 8'sd1;
                aEz_partial = acc_e - THRESHOLD;
            end else if (acc_e <= -THRESHOLD) begin
                Ez_partial  = -8'sd1;
                aEz_partial = acc_e + THRESHOLD;
            end else begin
                Ez_partial  = Ez_cur;
                aEz_partial = acc_e;
            end

            state <= S_E_SOURCE;
        end

        S_E_SOURCE: begin
            // Matches C exactly: the source cell is quantised a
            // SECOND time this timestep, after the source term is
            // added on top of whatever the normal update just did.
            if ((x == SX) && (y == SY)) begin
                acc_e2 = aEz_partial + sine_value;

                if (acc_e2 >= THRESHOLD) begin
                    Ez_final  = 8'sd1;
                    aEz_final = acc_e2 - THRESHOLD;
                end else if (acc_e2 <= -THRESHOLD) begin
                    Ez_final  = -8'sd1;
                    aEz_final = acc_e2 + THRESHOLD;
                end else begin
                    Ez_final  = Ez_partial;
                    aEz_final = acc_e2;
                end
            end else begin
                Ez_final  = Ez_partial;
                aEz_final = aEz_partial;
            end

            state <= S_E_WRITE;
        end

        S_E_WRITE: begin
            write_addr    <= addr_cur;
            write_data_ez <= Ez_final;
            acc_data_ez   <= aEz_final;
            write_to_h    <= 1'b0;
            write_enable  <= 1'b1;
            state         <= S_NEXT;
        end

        //=============================================================
        // NEXT CELL - note H and E phases sweep DIFFERENT bounds,
        // matching the C loops exactly (H starts at 0, E starts at 1).
        //=============================================================
        S_NEXT: begin
            write_enable <= 1'b0;

            if (h_phase) begin
                if (x == NX-2) begin
                    if (y == NY-2) begin
                        x       <= 7'd1;
                        y       <= 7'd1;
                        h_phase <= 1'b0;
                        state   <= S_E_READ_H;
                    end else begin
                        x <= 7'd0;
                        y <= y + 7'd1;
                        state <= S_H_READ;
                    end
                end else begin
                    x <= x + 7'd1;
                    state <= S_H_READ;
                end
            end else begin
                if (x == NX-2) begin
                    if (y == NY-2) begin
                        state <= S_WAIT;
                    end else begin
                        x <= 7'd1;
                        y <= y + 7'd1;
                        state <= S_E_READ_H;
                    end
                end else begin
                    x <= x + 7'd1;
                    state <= S_E_READ_H;
                end
            end
        end

        S_WAIT: begin
            state <= S_IDLE;
        end

        default: state <= S_IDLE;

        endcase
    end
end

//=====================================================================
// LED MONITOR
//
// Fields only ever take values -1, 0, +1 now, so the old
// abs_centre[15:14] magnitude-bit display no longer means anything
// (it would never differ from 00). Replaced with a direct sign
// indicator instead.
//=====================================================================

reg [12:0] centre_addr;
reg signed [7:0] centre_val;

always @(posedge clk) begin
    centre_addr <= SY * NX + SX;
    centre_val  <= Ez_mem[centre_addr];
end

assign LED[0] = (state != S_IDLE) && (state != S_INIT_CLEAR);
assign LED[1] = h_phase;
assign LED[2] = ~h_phase;
assign LED[3] = (source_phase != 16'd0);
assign LED[4] = (centre_val < 0);
assign LED[5] = (centre_val > 0);

endmodule