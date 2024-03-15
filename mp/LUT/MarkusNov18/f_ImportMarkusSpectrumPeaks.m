function [Current,lambda] = f_ImportMarkusSpectrumPeaks(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [CURRENT,MODE1,MODE2,MODE3] = IMPORTFILE(FILENAME) Reads data from text
%   file FILENAME for the default selection.
%
%   [CURRENT,MODE1,MODE2,MODE3] = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [Current,Mode1,Mode2,Mode3] = importfile('OSI_Peaks_4um_Temp20_132233.txt',1, 12);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/09 08:00:51

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));


for col=1:size(dataArray,2)
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
%'veer', keyboard

%% Replace non-numeric cells with NaN
%R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
%raw(R) = {NaN}; % Replace non-numeric cells


%% Allocate imported array to column variable names
%Current = cell2mat(raw(:, 1));
%for k=1:size(raw,2)
%lambda(k,:) = cell2mat(raw(:, k+1));
%end

Current = numericData(:, 1);
for k=1:size(raw,2)-1
lambdad(k,:) = numericData(:, k+1);
end

LA=lambdad';
for k=1:size(LA,1)
 fi=find(isnan(LA(k,:))==0);
 LA(k,fi)=fliplr(sort(LA(k,fi)));
end 

%'veer', keyboard
for k=1:size(LA,1)
 fi=find(isnan(LA(k,:))==0);
 lai=LA(k,fi);
 dl=-diff([0 lai]);
 LA(k,fi)=fliplr(sort(LA(k,fi)));
 if length(dl)>1
  fis=find(abs(dl)>.5*max(dl));
  if length(fis)>=1
   fibuo=fis;
   lad= LA(k,fi);
   LA(k,fi)=NaN*ones(size(fi));
   LA(k,1:length(fibuo))=lad(fibuo);
  end
 end
end 

lambda=LA;