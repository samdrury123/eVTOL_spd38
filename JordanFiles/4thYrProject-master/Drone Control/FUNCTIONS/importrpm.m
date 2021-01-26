function [time, RPM, thr] = importrpm(inputfile)
%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/jordaneriksen/Documents/Uni/Part 2B/Project/Drone Control/STATIC_TEST.txt
%
% Auto-generated by MATLAB on 30-Jan-2020 16:34:40

switch inputfile
    case 'RPM'
        filepath = "/Users/jordaneriksen/Documents/Uni/Part 2B/Project/Drone Control/STATIC_RPM.txt";
    case 'TEST'
        filepath = "/Users/jordaneriksen/Documents/Uni/Part 2B/Project/Drone Control/STATIC_TEST.txt";
    case 'THRUST'
        filepath = "/Users/jordaneriksen/Documents/Uni/Part 2B/Project/Drone Control/STATIC_THRUST.txt";
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["time", "RPM", "thr"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filepath, opts);

%% Convert to output type
time = tbl.time;
RPM = tbl.RPM;
thr = tbl.thr;

%% Clear temporary variables
clear opts tbl