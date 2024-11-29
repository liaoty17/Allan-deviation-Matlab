function [ProcessedData] = Data_Preprocessing(LiveData, IfRemoveDrift, IfRemoveOutliers, DataType)
% Data preprocessing. 

    % Copyright (C) 2024  Tangyin Liao
    % 
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <https://www.gnu.org/licenses/>.
    % 

% LiveData : the data using for calculate the allan deviation, should be 1-column.
% IfRemoveDrift : Whether to remove drift from data. Type:logical, 1 or 0.
% IfRemoveOutliers : Whether to remove outliers from data. Type:logical, 1 or 0.
% Datatype : input data type. Parameter type, "String"
    % "Phase" : Phase data.
    % "Frequency" : <default> Frequency data. 
% Output
    % ProcessedData: Data after processing.
% Update:2024-11-29, release: 1.0 Tangyin Liao
    % Add RemoveOutliers
    % Add RemoveDrift
    
%% Input parameters
    [rows, cols] = size(LiveData);
    if min(cols, rows) > 1
        error('Matlab::Arrays have multiple dimensions. Must be one dimension!')
    end
    if nargin < 1 || isempty(LiveData)
        error('Matlab::Data is empty or no parameter contribute to function!')
    elseif nargin == 1
        DataType = "Frequency";
        IfRemoveDrift = true;
        IfRemoveOutliers = true;
    elseif nargin == 2
        IfRemoveOutliers = true;
        DataType = "Frequency";
    elseif nargin ==3
        DataType = "Frequency";
    end
    
    %% main code
    N = length(LiveData);
    window = N/100;
    if IfRemoveDrift
        outData_nodrift = RemoveDrift(LiveData, DataType);
        if IfRemoveOutliers
            ProcessedData = RemoveOutliers(outData_nodrift, window);
        else
            ProcessedData = outData_nodrift;
        end
    end
end

function data_out = RemoveDrift(inputdata,DataType)
% Remove drift function, used by Data_Preprocessing
% Input:
    % inputdata: 1 column data. Phase or frequency time series.
    % DataType : "Frequency" or "Phase", use the default parameter "Frequency" if without input.
    if nargin == 1 
        DataType = "Frequency"; % Remove linear drift. It can be used in most cases.
    end
    if isempty(inputdata)
        error('Matlab::Data is empty or no parameter contribute to function!')
    end
    [rows, cols] = size(inputdata);
    if rows < cols
        inputdata = inputdata';
    end

    if strcmp(DataType, "Frequency")
        x = (1:1:length(inputdata))';
        kxb = polyfit(x, inputdata, 1);  % Linear drift.
        data_out = inputdata - polyval(kxb,x); 
    elseif strcmp(DataType, "Phase")
        x = (1:1:length(inputdata))';
        kxb = polyfit(x, inputdata, 2); % Remove second-order drift for phase data.
        data_out = inputdata - polyval(kxb,x);
    end
end

function out_data = RemoveOutliers(inputdata, window)
% Remove outliers, used by Data_Preprocessing
% Input:
    % inputdata: 1 column data. Phase or frequency or other time series.
    % window: move window , see <help rmoutliers>.
    if isempty(inputdata)
        error('Matlab::Data is empty or no parameter contribute to function!')
    end
    if nargin<2
        window = 20; % Move window length.
    end
    [rows, cols] = size(inputdata);
    if rows < cols
        inputdata = inputdata';
    end
    [out_data, TFrm] = rmoutliers(inputdata,"movmedian", window);
end