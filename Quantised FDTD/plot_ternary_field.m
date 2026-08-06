clear;
close all;
clc;

NX = 80;
NY = 80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select four snapshots here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steps = [180 260 340 480];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FONT_SIZE  = 14;
TITLE_SIZE = 16;

% FP sign threshold
SIGN_THRESHOLD = 0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read snapshot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = read_snapshot(filename,NX,NY)

    fid = fopen(filename,"r");

    if(fid < 0)
        error("Cannot open file: %s",filename);
    endif

    % Skip header
    fgetl(fid);
    fgetl(fid);

    A = fscanf(fid,"%f",[NY NX])';

    fclose(fid);

endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(100);
clf;

for n = 1:4

    step = steps(n);

    ternary_file = sprintf("snapshot_%04d.txt",step);
    fp_file      = sprintf("fp_snapshot_%04d.txt",step);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load fields
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    T = read_snapshot(ternary_file,NX,NY);
    F = read_snapshot(fp_file,NX,NY); F = sign(F);

    % Normalise floating point field
    F ./= max(abs(F(:)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spatial polarity agreement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Allow one-cell neighbourhood
    kernel = ones(3,3);

    FP_positive = conv2(F>0,kernel,"same") > 0;
    FP_negative = conv2(F<0,kernel,"same") > 0;

    pos = find(T == 1);
    neg = find(T == -1);

    matches = 0;
    total = length(pos)+length(neg);

    if(length(pos)>0)
        matches += sum(FP_positive(pos));
    endif

    if(length(neg)>0)
        matches += sum(FP_negative(neg));
    endif

    overlap = 100*matches/total;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot overlay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(2,2,n)

    imagesc(F)

    axis image
    set(gca,"YDir","normal")
    set(gca,"FontSize",FONT_SIZE)


    colormap(gray)

    hold on


    % Positive ternary field (+1)

    [row,col] = find(T == 1);

    plot(col,row,...
         'r.',...
         'MarkerSize',8)



    % Negative ternary field (-1)

    [row,col] = find(T == -1);

    plot(col,row,...
         'b.',...
         'MarkerSize',8)


    hold off



    title(sprintf("Step %d  (agreement %.1f%%)",...
          step,overlap),...
          "FontSize",TITLE_SIZE)

endfor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common colour scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%h = colorbar;

%set(h,"FontSize",FONT_SIZE)
