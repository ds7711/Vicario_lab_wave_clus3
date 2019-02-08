% Find the valid data files in the specified folder and store their path
% and name

function [FileName] = File_Detection(keyword, Starting_Path)
%% Select the folder storing the data and extract the filename of effective files; 
% The filename extracted are stored in cell array 'FileName'

% Select the directory where the data is stored and then put the file list
% in 'FileList' 
FileList = dir(Starting_Path);

Num_Files = size(FileList,1); % # of files in the folder
% FileName=cell(Num_Files,1);

% The first loop determines the size of the cell "FileName"
j=1;  % index for the data files (ID and trig file are not included)
for i = 1:Num_Files  % check every file inside the folder
    tempName = FileList(i,1).name;  % get the name file for below use
    judge=sum(regexpi(tempName,keyword));  % check if the filename contain the required keyword
    if judge
        j=j+1;
    end
end

FileName=cell(j-1,1);  % There're j-1 units in the data folder
j=1;  % index for the data files (ID and trig file are not included)
for i=1:Num_Files  % check every file inside the folder
    tempName=FileList(i,1).name;  % get the name file for below use
    judge=sum(regexpi(tempName,keyword));  % check if the filename contain the required keyword
    if judge
        FileName{j,1}=FileList(i,1).name;  % only the filename that contains the keyword is extracted into the cell arry
        j=j+1;
    end
end

if isempty(FileName{1,1})
    error('No Valid Data Files are detected! Please check the detection keyword or directory path!!!')
else
    disp('File path and file names are successfully extracted! Please Continue! ')
end


end