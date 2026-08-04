clear;
close all;
clc;

NX = 80;
NY = 80;

ternary_file = "snapshot_0480.txt";
fp_file      = "fp_snapshot_0480.txt";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read snapshot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = read_snapshot(filename,NX,NY)

    fid = fopen(filename,"r");

    fgetl(fid);
    fgetl(fid);

    A = fscanf(fid,"%f",[NY NX])';

    fclose(fid);

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = read_snapshot(ternary_file,NX,NY);
F = read_snapshot(fp_file,NX,NY);

F ./= max(abs(F(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FONT_SIZE = 16;       
TITLE_SIZE = 18;      
LINE_WIDTH = 2;

figure(1);
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field snapshots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1)

imagesc(T)

axis image
set(gca,"YDir","normal")
set(gca,"FontSize",FONT_SIZE)

colormap(gca,[0 0 1;1 1 1;1 0 0])

caxis([-1 1])

colorbar

title("Ternary field","FontSize",TITLE_SIZE)


subplot(2,3,2)

imagesc(F)

axis image
set(gca,"YDir","normal")
set(gca,"FontSize",FONT_SIZE)

colormap(gca,jet)

caxis([-1 1])

colorbar

title("Floating-point field","FontSize",TITLE_SIZE)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FT = abs(fftshift(fft2(T)));
FF = abs(fftshift(fft2(F)));

subplot(2,3,4)

imagesc(FT+1)

axis image
set(gca,"YDir","normal")
set(gca,"FontSize",FONT_SIZE)

colorbar

title("Spectrum (ternary)","FontSize",TITLE_SIZE)


subplot(2,3,5)

imagesc(FF+1)

axis image
set(gca,"YDir","normal")
set(gca,"FontSize",FONT_SIZE)

colorbar

title("Spectrum (floating-point)","FontSize",TITLE_SIZE)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radial spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cx = floor(size(FT,2)/2)+1;
cy = floor(size(FT,1)/2)+1;

R = floor(min(cx,cy));

specT = zeros(R,1);
specF = zeros(R,1);
count = zeros(R,1);


for i=1:size(FT,1)

    for j=1:size(FT,2)

        r = round(sqrt((i-cy)^2+(j-cx)^2));

        if(r>=1 && r<=R)

            specT(r) += FT(i,j);
            specF(r) += FF(i,j);
            count(r) += 1;

        endif

    endfor

endfor


specT ./= count;
specF ./= count;


subplot(2,3,[3 6])


plot(10*log10(abs(specT)/max(abs(specT))),...
     "LineWidth",LINE_WIDTH)

hold on

plot(10*log10(abs(specF)/max(abs(specF))),...
     "LineWidth",LINE_WIDTH)


grid on

set(gca,"FontSize",FONT_SIZE)

xlabel("Spatial frequency","FontSize",FONT_SIZE)

ylabel("Average spectral magnitude","FontSize",FONT_SIZE)

legend("Ternary","Floating-point",...
       "FontSize",FONT_SIZE)

title("Radially averaged spatial spectrum",...
      "FontSize",TITLE_SIZE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,"PaperPositionMode","auto")

print("Fig1_ch13_fdtd.pdf","-dpdf","-bestfit")

print("Fig1_ch13_fdtd.png","-dpng","-r300")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectral correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1 = FT(:)-mean(FT(:));
v2 = FF(:)-mean(FF(:));

spectral_corr = (v1'*v2)/(norm(v1)*norm(v2));


fprintf("\n");
fprintf("Spectral correlation = %.4f\n",spectral_corr);
fprintf("\n");
