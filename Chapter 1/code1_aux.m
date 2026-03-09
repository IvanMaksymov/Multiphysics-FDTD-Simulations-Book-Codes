% Octave script to load and plot Ez fields from 2D FDTD output
clc; close all; clear all;

% Choose step to plot
step = 300; filename = sprintf('Ez_step_%04d.dat', step);

% Load and plot the 2D field
Ez = load(filename);
figure; imagesc(Ez');        % transpose to match C indexing
axis equal tight; colorbar; xlabel('X'); ylabel('Y');
title(sprintf('Ez field at step %d', step)); colormap(jet);

