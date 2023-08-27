clc 
clear
close all
%%
distances = [11.4 23.2; 22 33.9; 31.7 36.2; 37.5 43.5; 29.3 36; 31.2 36; 20.5 29.8; 24.8 37.5; 14.8 20.4; 16.6 21.4];

innTDG = (distances(:,2) - distances(:,1))/340 * 1000;