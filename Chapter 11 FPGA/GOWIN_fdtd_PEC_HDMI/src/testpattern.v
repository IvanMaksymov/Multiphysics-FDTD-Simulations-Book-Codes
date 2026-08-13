// ---------------------------------------------------------------------
// Tile Display Test Pattern
// Colour format      : {B,G,R}
// ---------------------------------------------------------------------

module testpattern
(
    input              I_pxl_clk,
    input              I_rst_n,
    input      [2:0]   I_mode,

    input      [7:0]   I_single_r,
    input      [7:0]   I_single_g,
    input      [7:0]   I_single_b,

    input      [11:0]  I_h_total,
    input      [11:0]  I_h_sync,
    input      [11:0]  I_h_bporch,
    input      [11:0]  I_h_res,

    input      [11:0]  I_v_total,
    input      [11:0]  I_v_sync,
    input      [11:0]  I_v_bporch,
    input      [11:0]  I_v_res,

    input              I_hs_pol,
    input              I_vs_pol,

    input signed [15:0] I_fdtd_ez,

    output             O_de,
    output reg         O_hs,
    output reg         O_vs,

    output     [7:0]   O_data_r,
    output     [7:0]   O_data_g,
    output     [7:0]   O_data_b,

    output [11:0] O_fdtd_addr
);

//=====================================================================
// Colour definitions
// Format: {B,G,R}
//=====================================================================

localparam [23:0] WHITE   = {8'd255, 8'd255, 8'd255};
localparam [23:0] YELLOW  = {8'd0,   8'd255, 8'd255};
localparam [23:0] CYAN    = {8'd255, 8'd255, 8'd0};
localparam [23:0] GREEN   = {8'd0,   8'd255, 8'd0};
localparam [23:0] MAGENTA = {8'd255, 8'd0,   8'd255};
localparam [23:0] RED     = {8'd0,   8'd0,   8'd255};
localparam [23:0] BLUE    = {8'd255, 8'd0,   8'd0};
localparam [23:0] BLACK   = {8'd0,   8'd0,   8'd0};

//=====================================================================
// Display timing
//=====================================================================

localparam N = 5;

reg  [11:0] V_cnt;
reg  [11:0] H_cnt;

wire Pout_de_w;
wire Pout_hs_w;
wire Pout_vs_w;

reg [N-1:0] Pout_de_dn;
reg [N-1:0] Pout_hs_dn;
reg [N-1:0] Pout_vs_dn;

//=====================================================================
// Active-video counters
//=====================================================================

wire De_pos;
wire De_neg;
wire Vs_pos;

reg [11:0] De_vcnt;
reg [11:0] De_hcnt;

//=====================================================================
// Generate HS, VS and DE
//=====================================================================

always @(posedge I_pxl_clk or negedge I_rst_n)
begin
    if (!I_rst_n)
        V_cnt <= 12'd0;
    else begin
        if ((V_cnt >= (I_v_total - 1'b1)) &&
            (H_cnt >= (I_h_total - 1'b1)))
            V_cnt <= 12'd0;
        else if (H_cnt >= (I_h_total - 1'b1))
            V_cnt <= V_cnt + 1'b1;
        else
            V_cnt <= V_cnt;
    end
end

//---------------------------------------------------------------------

always @(posedge I_pxl_clk or negedge I_rst_n)
begin
    if (!I_rst_n)
        H_cnt <= 12'd0;
    else if (H_cnt >= (I_h_total - 1'b1))
        H_cnt <= 12'd0;
    else
        H_cnt <= H_cnt + 1'b1;
end

//---------------------------------------------------------------------

assign Pout_de_w =
       ((H_cnt >= (I_h_sync + I_h_bporch)) &&
        (H_cnt <= (I_h_sync + I_h_bporch + I_h_res - 1'b1))) &&

       ((V_cnt >= (I_v_sync + I_v_bporch)) &&
        (V_cnt <= (I_v_sync + I_v_bporch + I_v_res - 1'b1)));

//---------------------------------------------------------------------

assign Pout_hs_w =
       ~((H_cnt >= 12'd0) &&
         (H_cnt <= (I_h_sync - 1'b1)));

//---------------------------------------------------------------------

assign Pout_vs_w =
       ~((V_cnt >= 12'd0) &&
         (V_cnt <= (I_v_sync - 1'b1)));

//=====================================================================
// Pipeline DE / HS / VS
//=====================================================================

always @(posedge I_pxl_clk or negedge I_rst_n)
begin
    if (!I_rst_n) begin
        Pout_de_dn <= {N{1'b0}};
        Pout_hs_dn <= {N{1'b1}};
        Pout_vs_dn <= {N{1'b1}};
    end else begin
        Pout_de_dn <= {Pout_de_dn[N-2:0], Pout_de_w};
        Pout_hs_dn <= {Pout_hs_dn[N-2:0], Pout_hs_w};
        Pout_vs_dn <= {Pout_vs_dn[N-2:0], Pout_vs_w};
    end
end

assign O_de = Pout_de_dn[4];

//=====================================================================
// Output HS / VS
//=====================================================================

always @(posedge I_pxl_clk or negedge I_rst_n)
begin
    if (!I_rst_n) begin
        O_hs <= 1'b1;
        O_vs <= 1'b1;
    end else begin
        O_hs <= I_hs_pol ? ~Pout_hs_dn[3] : Pout_hs_dn[3];
        O_vs <= I_vs_pol ? ~Pout_vs_dn[3] : Pout_vs_dn[3];
    end
end

//=====================================================================
// Detect active-video boundaries
//=====================================================================

assign De_pos = !Pout_de_dn[1] & Pout_de_dn[0];
assign De_neg = Pout_de_dn[1] & !Pout_de_dn[0];
assign Vs_pos = !Pout_vs_dn[1] & Pout_vs_dn[0];

//=====================================================================
// Horizontal active-pixel counter
//=====================================================================

always @(posedge I_pxl_clk or negedge I_rst_n)
begin
    if (!I_rst_n)
        De_hcnt <= 12'd0;
    else if (De_pos)
        De_hcnt <= 12'd0;
    else if (Pout_de_dn[1])
        De_hcnt <= De_hcnt + 1'b1;
    else
        De_hcnt <= De_hcnt;
end

//=====================================================================
// Vertical active-line counter
//=====================================================================

always @(posedge I_pxl_clk or negedge I_rst_n)
begin
    if (!I_rst_n)
        De_vcnt <= 12'd0;
    else if (Vs_pos)
        De_vcnt <= 12'd0;
    else if (De_neg)
        De_vcnt <= De_vcnt + 1'b1;
    else
        De_vcnt <= De_vcnt;
end

//=====================================================================
// FDTD TILE ADDRESS 
// 64 x 64 grid on 1280x720 display
//=====================================================================

// Each tile = 20 pixels wide (1280/64 = 20)
// Each tile = 11 pixels tall (720/64 = 11.25, use integer 11)

wire [5:0] tile_x;
wire [5:0] tile_y;

// Simple integer division
assign tile_x = De_hcnt / 12'd20;
assign tile_y = De_vcnt / 12'd11;

// Clamp to 0-63 (safety)
wire [5:0] tile_x_clamped = (tile_x > 6'd63) ? 6'd63 : tile_x;
wire [5:0] tile_y_clamped = (tile_y > 6'd63) ? 6'd63 : tile_y;

// FDTD address = y * 64 + x
assign O_fdtd_addr = tile_y_clamped * 12'd64 + tile_x_clamped;

//=====================================================================
// FDTD FIELD COLOUR MAPPING
//=====================================================================

reg  [23:0] Tile_color;
wire [15:0] abs_value;
wire [7:0]  intensity;

localparam signed [15:0] MAX_VALUE = 16'sd200;

assign abs_value = I_fdtd_ez[15] ? -I_fdtd_ez : I_fdtd_ez;
assign intensity = (abs_value >= MAX_VALUE) ? 8'd255 : (abs_value * 8'd255) / MAX_VALUE;

always @* begin
    if (I_fdtd_ez > 16'sd0) begin
        // Positive: White to Red
        Tile_color = {8'd255, 8'd255 - intensity, 8'd255 - intensity};
    end else if (I_fdtd_ez < 16'sd0) begin
        // Negative: White to Blue
        Tile_color = {8'd255 - intensity, 8'd255 - intensity, 8'd255};
    end else begin
        // Zero: White
        Tile_color = 24'hFFFFFF;
    end
end

//=====================================================================
// RGB OUTPUT - Format {B,G,R}
//=====================================================================

assign O_data_r = Tile_color[7:0];
assign O_data_g = Tile_color[15:8];
assign O_data_b = Tile_color[23:16];

endmodule