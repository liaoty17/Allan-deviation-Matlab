function [Out_tau, TADEV_Output, Out_ErrBar] =Allan_TADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)
% Modified allan deviation calculation. 
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
% Update:2024-11-27 21:31, release: 2.0 Tangyin Liao
    % Code optimization to solve rounding errors.

%% Input parameters
    switch nargin
        case 1
            error('Matlab::Not enough input parameter arguments!')
        case 2
            Phase_or_Fre = "Frequency";
            len_m = floor(log2(length(LiveData)/2)); % The max length can be calculate.
            m = 2.^(0:1:len_m-1)';                   % Octave: Calculate using the exponential sequence of 2. Subtracting 1 is used to calculate the EDF.
            tau = tau0.*m;            
            Conf_interval = 0.683;
            Noise_type = "Auto";
        case 3          
            len_m = floor(log2(length(LiveData)/2)); 
            m = 2.^(0:1:len_m-1)';                     
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
    [rows, cols] = size(LiveData);
    if min(cols, rows) > 1
        error('Matlab::Arrays have multiple dimensions. Must be one dimension!')
    end
    if max(tau) > tau0*floor(length(LiveData)/2)
        error('Matlab::The maximum tau has exceeded the allowable limit.')
    end
    if length(LiveData) < 32
        error("Not enough data to calculate, should >32 points.")
    end
    %% Main code
    % Change array direction
    if cols > rows
        LiveData = LiveData'; % The reason for doing this is that when m=1, the indexing is prone to errors in N*1 array.
    end
    if strcmp(Phase_or_Fre, "Frequency")    
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        TADEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1);
        for i=1:1:len_tau
            N = length(LiveData);  
            LiveData_RUse_OS_Average = nan(N-m(i)+1, 1); % N-m+1
            LiveData_RUse_OS_Average(1) = mean(LiveData(1:m(i)));
            for ii = 2:1:N-m(i)+1
                LiveData_RUse_OS_Average(ii) = LiveData_RUse_OS_Average(ii-1) + (LiveData(ii+m(i)-1)-LiveData(ii-1))/m(i);
            end
            Delta_yData_bar = LiveData_RUse_OS_Average(1+m(i):N-m(i)+1)-LiveData_RUse_OS_Average(1:N-2*m(i)+1); % N-2m+1
            Delta_yData_bar_mod_mean = nan(N-3*m(i)+2, 1); % N-3m+2
            Delta_yData_bar_mod_mean(1) = mean(Delta_yData_bar(1:m(i)));
            for kk = 2:1:N-3*m(i)+2  
                Delta_yData_bar_mod_mean(kk) = Delta_yData_bar_mod_mean(kk-1) + (Delta_yData_bar(kk+m(i)-1)-Delta_yData_bar(kk-1))/m(i);
            end
            TADEV_Output(i,1) = sqrt(0.5*(mean(Delta_yData_bar_mod_mean.^2)))*Out_tau(i)/sqrt(3);
    
            % Calculate the error bar. 
            switch Noise_type
                case "Auto"
                    if floor(N/m(i)) >= 32
                        log_slope_alpha_approximate = Lag1_NoiseID(LiveData_RUse_OS_Average, 0, 2, Phase_or_Fre, m(i));
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
            edf_msig = EDF_Numerical_sim(log_slope_alpha_approximate,2,m(i),1, m(i), N+1);
            Out_ErrBar(i, 1:2) = TADEV_Output(i,1).*[sqrt(edf_msig/chi2inv(0.5+Conf_interval/2, edf_msig)), sqrt(edf_msig/chi2inv(0.5-Conf_interval/2, edf_msig))];
        end
            
    elseif strcmp(Phase_or_Fre, "Phase")
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        TADEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1);
        for j=1:1:len_tau
            N = length(LiveData);  % The phase data length.
            LiveData_RU_Diff = LiveData(1+m(j):N) - LiveData(1:N-m(j)); % N-m
            LiveData_RU_Diff2 = LiveData_RU_Diff(1+m(j):N-m(j)) - LiveData_RU_Diff(1:N-2*m(j)); % N-2m
            Delta_xData_mod_mean = nan(N-3*m(j)+1, 1); % N-3m+1
            Delta_xData_mod_mean(1) = mean(LiveData_RU_Diff2(1:m(j)));
            for kk =2:1:N-3*m(j)+1
                Delta_xData_mod_mean(kk) = Delta_xData_mod_mean(kk-1) + (LiveData_RU_Diff2(kk+m(j)-1)-LiveData_RU_Diff2(kk-1))/m(j);
            end         
            TADEV_Output(j,1) = sqrt(0.5*(mean(Delta_xData_mod_mean.^2)/((Out_tau(j))^2)))*Out_tau(j)/sqrt(3);

            % Calculate the error bar.    
            switch Noise_type
                case "Auto"
                    if floor(N/m(j)) >= 32
                        log_slope_alpha_approximate = Lag1_NoiseID(LiveData(1:N-m(j)), 1, 2, Phase_or_Fre, m(j)); % dmin should be 1 for phase data.
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
            edf_msig = EDF_Numerical_sim(log_slope_alpha_approximate,2, m(j), 1, m(j), N);
            Out_ErrBar(j, 1:2) = TADEV_Output(j,1).*[sqrt(edf_msig/chi2inv(0.5+Conf_interval/2, edf_msig)), sqrt(edf_msig/chi2inv(0.5-Conf_interval/2, edf_msig))];
        end
    else
        error("Matlab::The data type name is not permitted, 'Frequency' or 'Phase' instead.")
    end
end

