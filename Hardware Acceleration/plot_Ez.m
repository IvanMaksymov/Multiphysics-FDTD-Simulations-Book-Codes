% plot_Ez.m
%
% Load and plot the Ez field saved by the CUDA FDTD code.
% The file Ez.dat must be in the same directory or provide full path.

clear all;
close all;
clc;

% ------------------------------------------------------------
% Load Ez field
% ------------------------------------------------------------

fprintf("Loading Ez.dat...\n");

Ez = dlmread("Ez.dat");   % reads numeric matrix from text file

% Check size
[nrows, ncols] = size(Ez);
fprintf("Loaded field of size %d x %d\n", nrows, ncols);

% ------------------------------------------------------------
% Plot Ez as an image
% ------------------------------------------------------------

figure(1);
imagesc(Ez);              % plot as 2D image
axis equal tight;         % square pixels, no borders
colormap(jet);            % nice color map
colorbar;                 % show scale

title("Ez field");
xlabel("x index");
ylabel("y index");

% ------------------------------------------------------------
% Optional: plot a cross-section through the center
% ------------------------------------------------------------

figure(2);
mid = round(nrows/2);
plot(Ez(mid, :), "LineWidth", 2);
grid on;
title("Ez cross-section at mid y");
xlabel("x index");
ylabel("Ez value");

fprintf("Done.\n");

