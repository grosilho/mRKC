% Plot solutions 
clear;
clc;
close all;

folder = '../dist/Release/GNU_Version_10-MacOSX/Tests/ManyDiffusionTerms/';
sol_name = 'Platen';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

figure;
plot(t,y);