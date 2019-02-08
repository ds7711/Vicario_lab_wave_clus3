% check errors in recording filenames

%% clear all
clear all; close all; clc;

%% select the files to be processed
% old files probably cannot be batch processed due to lack of embedded
% duration files
% use File_Detection to get list of files to be processed

keyword = '.*smr';  % '.' any character; '*' any number of character
starting_path = 'D:\Google_Drive\Lab_Data\';
original_data_folder = uigetdir(starting_path, 'Select the Folder where the data is stored'); 
original_data_folder = [original_data_folder, '\']; % add a "\" for completeness
% new_folder = 'matrix_data\';
% new_path = [original_data_folder, new_folder];
% mkdir(new_path); % create a new folder to store the data

[filenames] = File_Detection(keyword, original_data_folder);
[fn_pairs, unpaired_files] = find_data_pairs(filenames);