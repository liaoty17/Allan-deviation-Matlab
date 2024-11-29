function [alpha, d] = Lag1_NoiseID(DataVector, dmin, dmax, DataType, m)
% Lag1 Autocorrelation noise type identification
 % 
    % Copyright (C) 2024 Tangyin Liao
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

% Input:
    % DataVector :  Imput vector of phase or frequency data.
    % dmin : the minimum order of differencing(default = 0).
    % dmax : the maximum order of differencing.
    % m : Overlapping factor.
% Output:
    % alpha: an estimate of the α of the dominant power law noise type
    % d: (optionally) the value of d

% Reference: W. J. Riley and C. A. Greenhall, "Power law noise identification using the lag 1 autocorrelation," 2004 18th (EFTF 2004), Guildford, 2004, pp. 576-580.
% Lag 1 autocorrelation. Edit by Tangyi Liao
% Update:2024-11-26,16:11 release: 2.0 Tangyin Liao
    % Update for overlapping samplling process.
   
    if nargin < 1 || isempty(DataVector)
        error('Matlab::Data is empty or no parameter contribute to function!')
    elseif nargin == 1
        dmin = 0;
        dmax = 2;
        DataType = "Frequency";
    elseif nargin == 2
        dmax = 2;
        DataType = "Frequency";
    elseif nargin == 3
        DataType = "Frequency";
    elseif nargin == 4
        m = 1; % default, no overlapped sampling.
    elseif nargin > 5
        error('Matlab::The number of input parameters exceeds the maximum allowable limit.')
    end
    [rows, cols] = size(DataVector);
    if cols > rows
        DataVector = DataVector'; 
    end
    %% main code
    Done = true;
    d = 0;         % Order of differencing
    tabele_Lag1_r1 = [0,   -1/2,  -1/2,  -2/3, -2/3, -3/4;...
                    0.7, -1/3,  -1/3,  -3/5, -3/5, -5/7;...
                    1,    0,     0,    -1/2, -1/2, -2/3;...
                    1,    0.7,   0.7,  -1/3, -1/3, -3/5;...
                    1,    1,     1,     0,     0,  -1/2];
    alpha_table = [2,1,0,-1,-2];

    while Done
       N = length(DataVector);
       DataVector_aver = mean(DataVector); 
       AutoPara_1 = sum((DataVector(1:N-m)-DataVector_aver).*(DataVector(1+m:N)-DataVector_aver));
       AutoPara_2 = sum((DataVector-DataVector_aver).^2);
       r1 = AutoPara_1/AutoPara_2;
       delta_judge = r1/(1+r1);
       if d >= dmin && (delta_judge<0.25 || d >= dmax)
           Done = false;
           % p = -round(2*delta_judge,0)-2*d;
       else
           DataVector = DataVector(1+m:N) - DataVector(1:N-m);
           d = d+1;
       end
    end
    if strcmp(DataType, "Frequency")
        new_array = r1- tabele_Lag1_r1(:,2*(d+1));
        [row, ~]= find(abs(new_array) == min(abs(new_array)));
        alpha = alpha_table(1, row);
        % alpha = p;
    elseif strcmp(DataType, "Phase")
       new_array = r1- tabele_Lag1_r1(:,2*d+1);
        [row, ~]= find(abs(new_array) == min(abs(new_array)));
        alpha = alpha_table(1, row);
       % alpha = p+2;
    end

end