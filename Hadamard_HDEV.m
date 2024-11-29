function [Out_tau, HDEV_Output, Out_ErrBar, alpha_memory] = Hadamard_HDEV(LiveData, tau0, Phase_or_Fre, tau, Conf_interval, Noise_type)
% Hadamard deviation calculation.
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

% LiveData : the data using for calculate the allan deviation, should be 1-D.
% tau0 : time interval.
% tau : Average time sequence.
% Phase_or_Fre : input data type. Parameter type, "String"
    % "Phase" : Phase data.
    % "Frequency" : Frequency data.
% Conf_interval : Confidence interval, using for calculate the error bar.
% Noise_type : the noise type of your data, if in doubt, "Auto" is recommended. Parameter type, "String";
    % "RWFM" : Random walk frequency modulation, h_{-2}
    % "FFM" : Flicker frequency modulation, h_{-1}
    % "WFM" : White frequency modulation, h_{0}
    % "FPM" : Flicker phase modulation, h_{1}
    % "WPM" : White phase modulation, h_{2}
% Output
    % Out_tau: Average time sequence.
    % ADEV_Output: Allan deviation(No Overlapped)
    % Out_ErrBar: 
        % Column 1: Min sigma
        % Column 2: max sigma
% Reference: Riley, W. and Howe, D. (2008), Handbook of Frequency Stability Analysis, Special Publication (NIST SP)...
% Update:2024-11-28, release: 1.0
    % Code optimization to solve rounding errors.

%% Input parameters
    switch nargin
        case 1
            error('Matlab::Not enough input arguments!')
        case 2
            Phase_or_Fre = "Frequency";
            len_m = floor(log2(length(LiveData)/2)); % The max length can be calculate.
            m = 2.^(0:1:len_m-2)';                   % Octave: Calculate using the exponential sequence of 2. Subtracting 2 is used to calculate the EDF.
            tau = tau0.*m;            
            Conf_interval = 0.683;
            Noise_type = "Auto";
        case 3          
            len_m = floor(log2(length(LiveData)/2)); 
            m = 2.^(0:1:len_m-2)';                     
            tau = tau0.*m;
            Conf_interval = 0.683;
            Noise_type = "Auto";
        case 4       
            Conf_interval = 0.683;
            Noise_type = "Auto";
        case 5      
            Noise_type = "Auto";
    end
    if nargin < 1 || isempty(LiveData)
        error('Matlab::Data is empty or no parameter contribute to function!')
    end
    if min(size(LiveData)) > 1
        error('Matlab::Arrays have multiple dimensions. Must be one dimension!')
    end
    if max(tau) > tau0*floor(length(LiveData)/2)
        error('Matlab::The maximum tau has exceeded the allowable limit.')
    end    
    if length(LiveData) < 32
        error("Not enough data to calculate, should >32 points.")
    end
    %% Main code
    if strcmp(Phase_or_Fre, "Frequency")
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        HDEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1); 
        for i=1:1:len_tau
            M = floor(length(LiveData)/m(i));  % allan total blocks to average.
            LiveData_RealUse = LiveData(1:m(i)*M);
            LiveData_RealUse_Reshape = reshape(LiveData_RealUse, m(i),[]);
            LiveData_RUse_Re_Average = mean(LiveData_RealUse_Reshape, 1);
            N = length(LiveData_RUse_Re_Average);
            LiveData_Aver_Diff = LiveData_RUse_Re_Average(2:N)-LiveData_RUse_Re_Average(1:N-1);
            LiveData_Aver_Diff2 = LiveData_Aver_Diff(2:N-1) - LiveData_Aver_Diff(1:N-2);
            LiveData_Average_SQ = mean(LiveData_Aver_Diff2.^2);       
            HDEV_Output(i,1) = sqrt(LiveData_Average_SQ/6);
    
            % Calculate the error bar. 
            switch Noise_type
                case "Auto"
                    if M >= 32
                        log_slope_alpha_approximate = Lag1_NoiseID(LiveData_RUse_Re_Average, 0, 2, Phase_or_Fre, 1);
                        alpha_memory(i,1) = log_slope_alpha_approximate;
                    else
                        log_slope_alpha_approximate = alpha_memory(i-1,1);
                        alpha_memory(i,1) = log_slope_alpha_approximate;
                    end
                case "RWFM"
                    log_slope_alpha_approximate = -2;
                case "FFM"
                    log_slope_alpha_approximate = -1;
                case "WFM"
                    log_slope_alpha_approximate = 0;
                case "FPM"
                    log_slope_alpha_approximate = 1;
                case "WPM"
                    log_slope_alpha_approximate = 2;
            end            
            edf_hdev = EDF_Numerical_sim(log_slope_alpha_approximate,3,m(i),m(i),1,m(i)*M+1);
            Out_ErrBar(i, 1:2) = HDEV_Output(i,1).*[sqrt(edf_hdev/chi2inv(0.5+Conf_interval/2, edf_hdev)), sqrt(edf_hdev/chi2inv(0.5-Conf_interval/2, edf_hdev))];
        end
            
    elseif strcmp(Phase_or_Fre, "Phase")
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        HDEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1); 
        for j=1:1:len_tau
            N = floor((length(LiveData)-1)/m(j));      % The phase data length is subtracted by 1 to determine the number of average blocks.
            LiveData_Aver_Diff = LiveData(1:m(j)*N+1); % for phase data is add 1 to get the really used data number.
            LiveData_RealUse_Select = LiveData_Aver_Diff(1:m(j):end);
            LiveData_RealUse_Sel_Diff = LiveData_RealUse_Select(2:end)-LiveData_RealUse_Select(1:end-1); % 
            LiveData_Diff_Aver = LiveData_RealUse_Sel_Diff(2:end)-LiveData_RealUse_Sel_Diff(1:end-1); % 
            Delta_yData_bar = (LiveData_Diff_Aver(2:end)-LiveData_Diff_Aver(1:end-1)).^2;
            HDEV_Output(j,1) = sqrt((mean(Delta_yData_bar)/(Out_tau(j)^2)/6));

            % Calculate the error bar.    
            switch Noise_type
                case "Auto"
                    if length(LiveData_RealUse_Select) >= 32             
                        log_slope_alpha_approximate = Lag1_NoiseID(LiveData_RealUse_Select, 1, 2, Phase_or_Fre, 1); 
                        alpha_memory(j,1) = log_slope_alpha_approximate;
                    else
                        log_slope_alpha_approximate = alpha_memory(j-1,1);
                        alpha_memory(j,1) = log_slope_alpha_approximate;
                    end
                case "RWFM"
                    log_slope_alpha_approximate = -2;
                case "FFM"
                    log_slope_alpha_approximate = -1;
                case "WFM"
                    log_slope_alpha_approximate = 0;
                case "FPM"
                    log_slope_alpha_approximate = 1;
                case "WPM"
                    log_slope_alpha_approximate = 2;
            end            
            edf_hdev = EDF_Numerical_sim(log_slope_alpha_approximate,3,m(j),m(j),1,m(j)*N+1);
            Out_ErrBar(j, 1:2) = HDEV_Output(j,1).*[sqrt(edf_hdev/chi2inv(0.5+Conf_interval/2, edf_hdev)), sqrt(edf_hdev/chi2inv(0.5-Conf_interval/2, edf_hdev))];
        end
        
    else
        error("Matlab::The data type name is not permitted, 'Frequency' or 'Phase' instead.")
    end
end

