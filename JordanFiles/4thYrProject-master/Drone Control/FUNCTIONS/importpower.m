function [TIMESTAMPms, CURRENT, VOLTAGE] = importpower()
%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/jordaneriksen/Documents/Uni/Part 2B/Project/Drone Control/Logs/Current test1_compressed.txt
%
% Auto-generated by MATLAB on 03-Feb-2020 16:59:23

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["TIMESTAMPms", "CURRENT", "VOLTAGE"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("/Users/jordaneriksen/Documents/Uni/Part 2B/Project/Drone Control/Logs/power_compressed.txt", opts);

%% Convert to output type
TIMESTAMPms = tbl.TIMESTAMPms;
CURRENT = tbl.CURRENT;
VOLTAGE = tbl.VOLTAGE;

%% Clear temporary variables
clear opts tbl
end