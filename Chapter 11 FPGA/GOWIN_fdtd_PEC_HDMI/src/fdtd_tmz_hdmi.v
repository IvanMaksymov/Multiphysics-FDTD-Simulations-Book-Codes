`timescale 1ns / 1ps

module fdtd_tmz
(
    input  wire               clk,
    input  wire               rst,

    input  wire [11:0]        display_addr,
    output wire signed [15:0] display_ez,

    output wire [5:0]         LED
);

parameter NX   = 64;
parameter NY   = 64;
parameter SIZE = NX * NY;

//=====================================================================
// FDTD PARAMETERS
//=====================================================================

parameter STEP_DELAY = 26'd20000000;

// C = 0.25
localparam signed [15:0] C = 16'sd1024;

// Strong source
localparam signed [15:0] SOURCE_AMPLITUDE = 16'sd2000;
localparam [7:0] PHASE_INCREMENT = 8'd4;

//=====================================================================
// FIELD MEMORIES
//=====================================================================

(* ram_style = "block" *)
reg signed [15:0] Ez_mem [0:SIZE-1];

(* ram_style = "block" *)
reg signed [15:0] Hx_mem [0:SIZE-1];

(* ram_style = "block" *)
reg signed [15:0] Hy_mem [0:SIZE-1];

//=====================================================================
// DISPLAY READ
//=====================================================================

reg [11:0] display_addr_reg;
reg signed [15:0] display_ez_reg;

always @(posedge clk) begin
    display_addr_reg <= display_addr;
    display_ez_reg   <= Ez_mem[display_addr_reg];
end

assign display_ez = display_ez_reg;

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

reg [7:0] source_phase;
reg signed [15:0] sine_value;

always @(posedge clk) begin
    if (rst)
        source_phase <= 8'd0;
    else if (step)
        source_phase <= source_phase + PHASE_INCREMENT;
end

always @* begin
    case (source_phase[7:4])
        4'd0:  sine_value =  16'sd0;
        4'd1:  sine_value =  16'sd764;
        4'd2:  sine_value =  16'sd1414;
        4'd3:  sine_value =  16'sd1848;
        4'd4:  sine_value =  16'sd2000;
        4'd5:  sine_value =  16'sd1848;
        4'd6:  sine_value =  16'sd1414;
        4'd7:  sine_value =  16'sd764;
        4'd8:  sine_value =  16'sd0;
        4'd9:  sine_value = -16'sd764;
        4'd10: sine_value = -16'sd1414;
        4'd11: sine_value = -16'sd1848;
        4'd12: sine_value = -16'sd2000;
        4'd13: sine_value = -16'sd1848;
        4'd14: sine_value = -16'sd1414;
        4'd15: sine_value = -16'sd764;
        default: sine_value = 16'sd0;
    endcase
end

//=====================================================================
// FDTD CONTROLLER STATES
//=====================================================================

localparam [3:0]
    S_IDLE      = 4'd0,
    S_H_READ    = 4'd1,
    S_H_WAIT    = 4'd2,
    S_H_CALC    = 4'd3,
    S_H_WRITE   = 4'd4,
    S_E_READ_H  = 4'd5,
    S_E_WAIT_H  = 4'd6,
    S_E_READ_EZ = 4'd7,
    S_E_WAIT_EZ = 4'd8,
    S_E_CALC    = 4'd9,
    S_E_SOURCE  = 4'd10,
    S_E_WRITE   = 4'd11,
    S_NEXT      = 4'd12,
    S_WAIT      = 4'd13;

reg [3:0] state;
reg [5:0] x, y;
reg h_phase;

// Safer address generation
wire [5:0] x_m1 = (x == 6'd0)     ? 6'd0     : (x - 6'd1);
wire [5:0] x_p1 = (x == NX-1)     ? (NX-1)   : (x + 6'd1);
wire [5:0] y_m1 = (y == 6'd0)     ? 6'd0     : (y - 6'd1);
wire [5:0] y_p1 = (y == NY-1)     ? (NY-1)   : (y + 6'd1);

wire [11:0] addr_cur = y    * NX + x;
wire [11:0] addr_r   = y    * NX + x_p1;
wire [11:0] addr_u   = y_p1 * NX + x;
wire [11:0] addr_l   = y    * NX + x_m1;
wire [11:0] addr_d   = y_m1 * NX + x;

//=====================================================================
// REGISTERS FOR CALCULATIONS
//=====================================================================

reg signed [15:0] Ez_cur, Ez_r, Ez_u;
reg signed [15:0] Hx_cur, Hx_d;
reg signed [15:0] Hy_cur, Hy_l;
reg signed [15:0] Hx_new, Hy_new;
reg signed [15:0] Ez_new;
reg signed [31:0] curl;
reg signed [31:0] curl_x, curl_y;

//=====================================================================
// WRITE CONTROL - SEPARATE FROM READS
//=====================================================================

reg [11:0] write_addr;
reg signed [15:0] write_data_ez;
reg signed [15:0] write_data_hx;
reg signed [15:0] write_data_hy;
reg write_enable;
reg write_to_h;

//=====================================================================
// BRAM WRITE - SEPARATE PROCESS
//=====================================================================

always @(posedge clk) begin
    if (write_enable) begin
        if (write_to_h) begin
            Hx_mem[write_addr] <= write_data_hx;
            Hy_mem[write_addr] <= write_data_hy;
        end else begin
            Ez_mem[write_addr] <= write_data_ez;
        end
    end
end

//=====================================================================
// FDTD STATE MACHINE - FIXED FOR GOWIN BRAM
//=====================================================================

always @(posedge clk) begin
    if (rst) begin
        state <= S_IDLE;
        x <= 6'd1;
        y <= 6'd1;
        h_phase <= 1'b1;
        write_enable <= 1'b0;
    end else begin
        // Default: no write
        write_enable <= 1'b0;
        
        case (state)

        S_IDLE: begin
            if (step) begin
                x <= 6'd1;
                y <= 6'd1;
                h_phase <= 1'b1;
                state <= S_H_READ;
            end
        end

        //=============================================================
        // H FIELD UPDATE  (now with proper BRAM latency)
        //=============================================================
        S_H_READ: begin
            // Address is already stable (set in S_NEXT).
            // Just wait one cycle for the registered BRAM output.
            state <= S_H_WAIT;
        end

        S_H_WAIT: begin
            // Now the data is valid – capture it
            Ez_cur <= Ez_mem[addr_cur];
            Ez_r   <= Ez_mem[addr_r];
            Ez_u   <= Ez_mem[addr_u];
            Hx_cur <= Hx_mem[addr_cur];
            Hy_cur <= Hy_mem[addr_cur];
            state  <= S_H_CALC;
        end

        S_H_CALC: begin
            curl_x  = (Ez_u - Ez_cur) * C;
            curl_y  = (Ez_r - Ez_cur) * C;
            Hx_new <= Hx_cur - (curl_x >>> 12);
            Hy_new <= Hy_cur + (curl_y >>> 12);
            state  <= S_H_WRITE;
        end

        S_H_WRITE: begin
            write_addr    <= addr_cur;
            write_data_hx <= Hx_new;
            write_data_hy <= Hy_new;
            write_to_h    <= 1'b1;
            write_enable  <= 1'b1;
            state         <= S_NEXT;
        end

        //=============================================================
        // E FIELD UPDATE  (also with proper latency)
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
            Ez_cur <= Ez_mem[addr_cur];
            state  <= S_E_CALC;
        end

        S_E_CALC: begin
            curl   = (Hy_cur - Hy_l) - (Hx_cur - Hx_d);
            Ez_new <= Ez_cur + ((curl * C) >>> 12);
            state  <= S_E_SOURCE;
        end

        S_E_SOURCE: begin
            if ((x == 6'd32) && (y == 6'd32))
                Ez_new <= Ez_new + sine_value;
            state <= S_E_WRITE;
        end

        S_E_WRITE: begin
            write_addr   <= addr_cur;
            write_data_ez<= Ez_new;
            write_to_h   <= 1'b0;
            write_enable <= 1'b1;
            state        <= S_NEXT;
        end

        //=============================================================
        // NEXT CELL
        //=============================================================
        S_NEXT: begin
            // Clear write enable
            write_enable <= 1'b0;
            
            if (h_phase) begin
                // H sweep
                if (x == 6'd62) begin
                    if (y == 6'd62) begin
                        x <= 6'd1;
                        y <= 6'd1;
                        h_phase <= 1'b0;
                        state <= S_E_READ_H;
                    end else begin
                        x <= 6'd1;
                        y <= y + 6'd1;
                        state <= S_H_READ;
                    end
                end else begin
                    x <= x + 6'd1;
                    state <= S_H_READ;
                end
            end else begin
                // E sweep
                if (x == 6'd62) begin
                    if (y == 6'd62) begin
                        state <= S_WAIT;
                    end else begin
                        x <= 6'd1;
                        y <= y + 6'd1;
                        state <= S_E_READ_H;
                    end
                end else begin
                    x <= x + 6'd1;
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
//=====================================================================

reg [11:0] centre_addr;
reg signed [15:0] centre_val;

always @(posedge clk) begin
    centre_addr <= 32 * NX + 32;
    centre_val  <= Ez_mem[centre_addr];
end

wire signed [15:0] abs_centre = centre_val[15] ? -centre_val : centre_val;

assign LED[0] = (state != S_IDLE);
assign LED[1] = h_phase;
assign LED[2] = ~h_phase;
assign LED[3] = (source_phase != 8'd0);
assign LED[5:4] = abs_centre[15:14];

endmodule