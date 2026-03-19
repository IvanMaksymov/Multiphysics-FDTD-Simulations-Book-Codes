clc;clear all; close all;

data = load('field.dat');

for k = 1:size(data,1)
    plot(data(k,:));
    axis([0 2000 -1 1]);
    drawnow;
end
