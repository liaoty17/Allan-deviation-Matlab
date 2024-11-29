function [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = AllanTools_Matlab(DeviationType, LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)
% Allan deviation and Hadamard deviation package tools.

    % Copyright (C) <2024>  <Tangyin Liao>
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

% Input
% DeviationType: <String>
    % "ADEV" : standard Allan deviation.
    % "OADEV" : overlapping Allan deviation
    % "MADEV" : modified Allan deviation.
    % "TADEV" : time Allan deviation.
    % "HDEV" : Hadamard deviation.
    % "OHDEV" : overlapping Hadamard deviation.
    % "MHDEV" : modified Hadamard deviation.
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
    % alpha_memory: Estimation of noise type at the average time tau.

% V2.0 2024-11-29
    % 
    %% Input parameters
    switch nargin
        case 2
            error('Matlab::Not enough input arguments!')
        case 3
            Phase_or_Fre = "Frequency";
            len_m = floor(log2(length(LiveData)/2)); % The max length can be calculate.
            m = 2.^(0:1:len_m-2)';                   % Octave: Calculate using the exponential sequence of 2. Subtracting 2 is used to calculate the EDF.
            tau = tau0.*m;            
            Conf_interval = 0.683;
            Noise_type = "Auto";
        case 4          
            len_m = floor(log2(length(LiveData)/2)); 
            m = 2.^(0:1:len_m-2)';                     
            tau = tau0.*m;
            Conf_interval = 0.683;
            Noise_type = "Auto";
        case 5      
            Conf_interval = 0.683;
            Noise_type = "Auto";
        case 6      
            Noise_type = "Auto";
    end
    if nargin < 3 || isempty(LiveData)
        error('Matlab::Data is empty or not enough parameter contribute to function!')
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
    %% main code
    switch DeviationType
        case "ADEV"
            [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = Allan_ADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type);
        case "OADEV"
            [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = Allan_OADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type);
        case "MADEV"
            [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = Allan_MADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type);
        case "TADEV"
            [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = Allan_TADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type);
        case "HDEV"
            [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = Hadamard_HDEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type);
        case "OHDEV"
            [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = Hadamard_OHDEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type);
        case "MHDEV"
            [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = Hadamard_MHDEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type);
        otherwise
            disp("'ADEV' : standard Allan deviation.")
            disp("'OADEV' : overlapping Allan deviation")
            disp("'MADEV' : modified Allan deviation.")
            disp("'TADEV' : time Allan deviation.")
            disp("'HDEV' : Hadamard deviation.")
            disp("'OHDEV' : overlapping Hadamard deviation.")
            disp("'MHDEV' : modified Hadamard deviation.")
            error("Matlab::Chose the correct deviation model.")
      
    end
end


function [Out_tau, ADEV_Output, Out_ErrBar, alpha_memory] = Allan_ADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)
% Standard allan deviation calculation.Edit by Tangyi Liao
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
% Update:2024-11-27, release: 2.0
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
        ADEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1); 
        for i=1:1:len_tau
            M = floor(length(LiveData)/m(i));  % allan total blocks to average.
            LiveData_RealUse = LiveData(1:m(i)*M);
            LiveData_RealUse_Reshape = reshape(LiveData_RealUse, m(i),[]);
            LiveData_RUse_Re_Average = mean(LiveData_RealUse_Reshape, 1);
            Delta_yData_bar = (LiveData_RUse_Re_Average(2:end)-LiveData_RUse_Re_Average(1:end-1)).^2;
            ADEV_Output(i,1) = sqrt(0.5*(mean(Delta_yData_bar)));
    
            % Calculate the error bar. 
            switch Noise_type
                case "Auto"
                    if length(LiveData_RUse_Re_Average) >= 32
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
            edf_sig2 = EDF_Numerical_sim(log_slope_alpha_approximate,2,m(i),m(i),1,m(i)*M+1);
            Out_ErrBar(i, 1:2) = ADEV_Output(i,1).*[sqrt(edf_sig2/chi2inv(0.5+Conf_interval/2, edf_sig2)), sqrt(edf_sig2/chi2inv(0.5-Conf_interval/2, edf_sig2))];
        end

            
    elseif strcmp(Phase_or_Fre, "Phase")
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        ADEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1); 
        for j=1:1:len_tau
            N = floor((length(LiveData)-1)/m(j));  % The phase data length is subtracted by 1 to determine the number of average blocks.
            LiveData_RealUse = LiveData(1:m(j)*N+1); % for phase data is add 1 to get the really used data number.
            LiveData_RealUse_Select = LiveData_RealUse(1:m(j):end);
            LiveData_RealUse_Sel_Diff = LiveData_RealUse_Select(2:end)-LiveData_RealUse_Select(1:end-1);
            Delta_yData_bar = (LiveData_RealUse_Sel_Diff(2:end)-LiveData_RealUse_Sel_Diff(1:end-1)).^2;
            ADEV_Output(j,1) = sqrt(0.5*(mean(Delta_yData_bar)/((Out_tau(j))^2)));

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
            edf_sig2 = EDF_Numerical_sim(log_slope_alpha_approximate,2,m(j),m(j),1,m(j)*N+1);
            Out_ErrBar(j, 1:2) = ADEV_Output(j,1).*[sqrt(edf_sig2/chi2inv(0.5+Conf_interval/2, edf_sig2)), sqrt(edf_sig2/chi2inv(0.5-Conf_interval/2, edf_sig2))];
        end
        
    else
        error("Matlab::The data type name is not permitted, 'Frequency' or 'Phase' instead.")
    end
end

function [Out_tau, OADEV_Output, Out_ErrBar, alpha_memory] = Allan_OADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)
% Overlapping allan deviation calculation. Edit by Tangyi Liao
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
% Update:2024-11-27, release: 2.0
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
        LiveData = LiveData'; 
    end
    if strcmp(Phase_or_Fre, "Frequency")    
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        OADEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1);
        for i=1:1:len_tau
            N = length(LiveData);
            LiveData_RUse_OS_Average = nan(N-m(i)+1, 1);
            LiveData_RUse_OS_Average(1) = mean(LiveData(1:m(i)));
            for ii = 2:1:N-m(i)+1
                LiveData_RUse_OS_Average(ii) = LiveData_RUse_OS_Average(ii-1) + (LiveData(ii+m(i)-1)-LiveData(ii-1))/m(i);
            end
            Delta_yData_bar_Diff_SQ = (LiveData_RUse_OS_Average(1+m(i):N-m(i)+1)-LiveData_RUse_OS_Average(1:N-2*m(i)+1)).^2;
            OADEV_Output(i,1) = sqrt(0.5*(mean(Delta_yData_bar_Diff_SQ)));

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
            edf_osig = EDF_Numerical_sim(log_slope_alpha_approximate,2,m(i),m(i),m(i),N+1);
            Out_ErrBar(i, 1:2) = OADEV_Output(i,1).*[sqrt(edf_osig/chi2inv(0.5+Conf_interval/2, edf_osig)), sqrt(edf_osig/chi2inv(0.5-Conf_interval/2, edf_osig))];
        end
            
    elseif strcmp(Phase_or_Fre, "Phase")
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        OADEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1);
        for j=1:1:len_tau
            N = length(LiveData);  % The phase data length.
            LiveData_RU_Diff = LiveData(1+m(j):N) - LiveData(1:N-m(j)); % N-m
            LiveData_RU_Diff2_SQ = (LiveData_RU_Diff(1+m(j):N-m(j)) - LiveData_RU_Diff(1:N-2*m(j))).^2; % N-2m   
            OADEV_Output(j,1) = sqrt(0.5*(mean(LiveData_RU_Diff2_SQ)/((Out_tau(j))^2)));

            % Calculate the error bar.    
            switch Noise_type
                case "Auto"
                    if floor(N/m(j)) >= 32
                        log_slope_alpha_approximate = Lag1_NoiseID(LiveData, 1, 2, Phase_or_Fre, m(j)); % dmin should be 1.
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
            edf_osig = EDF_Numerical_sim(log_slope_alpha_approximate,2,m(j),m(j),m(j), N);
            Out_ErrBar(j, 1:2) = OADEV_Output(j,1).*[sqrt(edf_osig/chi2inv(0.5+Conf_interval/2, edf_osig)), sqrt(edf_osig/chi2inv(0.5-Conf_interval/2, edf_osig))];
        end
    else
        error("Matlab::The data type name is not permitted, 'Frequency' or 'Phase' instead.")
    end
end

function [Out_tau, MADEV_Output, Out_ErrBar, alpha_memory] = Allan_MADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)

% Modified allan deviation calculation. Edit by Tangyi Liao
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
% Update:2024-11-27 21:31, release: 2.0
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
        MADEV_Output  = nan(len_tau,1);
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
            MADEV_Output(i,1) = sqrt(0.5*(mean(Delta_yData_bar_mod_mean.^2)));
    
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
            edf_msig = EDF_Numerical_sim(log_slope_alpha_approximate,2, m(i), 1, m(i), N+1);
            Out_ErrBar(i, 1:2) = MADEV_Output(i,1).*[sqrt(edf_msig/chi2inv(0.5+Conf_interval/2, edf_msig)), sqrt(edf_msig/chi2inv(0.5-Conf_interval/2, edf_msig))];
        end
            
    elseif strcmp(Phase_or_Fre, "Phase")
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        MADEV_Output  = nan(len_tau,1);
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
            MADEV_Output(j,1) = sqrt(0.5*(mean(Delta_xData_mod_mean.^2)/((Out_tau(j))^2)));

            % Calculate the error bar.    
            switch Noise_type
                case "Auto"
                    if floor(N/m(j)) >= 32
                        log_slope_alpha_approximate = Lag1_NoiseID(LiveData, 1, 2, Phase_or_Fre, m(j)); % dmin should be 1 for phase data.
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
            Out_ErrBar(j, 1:2) = MADEV_Output(j,1).*[sqrt(edf_msig/chi2inv(0.5+Conf_interval/2, edf_msig)), sqrt(edf_msig/chi2inv(0.5-Conf_interval/2, edf_msig))];
        end
    else
        error("Matlab::The data type name is not permitted, 'Frequency' or 'Phase' instead.")
    end
end

function [Out_tau, TADEV_Output, Out_ErrBar] =Allan_TADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)
% Modified allan deviation calculation. Edit by Tangyi Liao
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
% Update:2024-11-27 21:31, release: 2.0
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

function [Out_tau, HDEV_Output, Out_ErrBar, alpha_memory] = Hadamard_HDEV(LiveData, tau0, Phase_or_Fre, tau, Conf_interval, Noise_type)
% Hadamard deviation calculation.Edit by Tangyi Liao
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

function [Out_tau, MHDEV_Output, Out_ErrBar, alpha_memory] = Hadamard_MHDEV(LiveData, tau0, Phase_or_Fre, tau, Conf_interval, Noise_type)
% Moidified Hadamard deviation calculation.Edit by Tangyi Liao
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
        MHDEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1); 
        for i=1:1:len_tau
            N = length(LiveData);  
            LiveData_RUse_OS_Average = nan(N-m(i)+1, 1); % N-m+1
            LiveData_RUse_OS_Average(1) = mean(LiveData(1:m(i)));
            for ii = 2:1:N-m(i)+1
                LiveData_RUse_OS_Average(ii) = LiveData_RUse_OS_Average(ii-1) + (LiveData(ii+m(i)-1)-LiveData(ii-1))/m(i);
            end
            Delta_yData_bar_Diff = LiveData_RUse_OS_Average(1+m(i):N-m(i)+1)-LiveData_RUse_OS_Average(1:N-2*m(i)+1); % N-2m+1
            Delta_yData_bar_Diff2 = Delta_yData_bar_Diff(1+m(i):N-2*m(i)+1)-Delta_yData_bar_Diff(1:N-3*m(i)+1);      % N-3m+1
            Delta_yData_bar_mod_mean = nan(N-4*m(i)+2, 1); % N-4m+2
            Delta_yData_bar_mod_mean(1) = mean(Delta_yData_bar_Diff2(1:m(i)));
            for kk = 2:1:N-4*m(i)+2  
                Delta_yData_bar_mod_mean(kk) = Delta_yData_bar_mod_mean(kk-1) + (Delta_yData_bar_Diff2(kk+m(i)-1)-Delta_yData_bar_Diff2(kk-1))/m(i);
            end       
            MHDEV_Output(i,1) = sqrt(mean(Delta_yData_bar_mod_mean.^2)/6);
    
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
            edf_mhdev = EDF_Numerical_sim(log_slope_alpha_approximate, 3, m(i), 1, m(i), N+1);
            Out_ErrBar(i, 1:2) = MHDEV_Output(i,1).*[sqrt(edf_mhdev/chi2inv(0.5+Conf_interval/2, edf_mhdev)), sqrt(edf_mhdev/chi2inv(0.5-Conf_interval/2, edf_mhdev))];
        end
            
    elseif strcmp(Phase_or_Fre, "Phase")
        len_tau = length(tau);
        m = floor(tau./tau0)';
        Out_tau = m.*tau0;
        MHDEV_Output  = nan(len_tau,1);
        Out_ErrBar = nan(len_tau,2); % Initial
        alpha_memory = nan(len_tau,1); 
        for j=1:1:len_tau
            N = length(LiveData);  % The phase data length.
            LiveData_to_yBar  = LiveData(1+m(j):N) - LiveData(1:N-m(j)); % N-m
            LiveData_yBar_Diff = LiveData_to_yBar(1+m(j):N-m(j)) - LiveData_to_yBar(1:N-2*m(j));        % N-2m
            LiveData_yBar_Diff2 = LiveData_yBar_Diff(1+m(j):N-2*m(j)) - LiveData_yBar_Diff(1:N-3*m(j)); % N-3m
            Delta_xData_mod_mean = nan(N-4*m(j)+1, 1); % N-4m+1
            Delta_xData_mod_mean(1) = mean(LiveData_yBar_Diff2(1:m(j)));
            for kk =2:1:N-4*m(j)+1
                Delta_xData_mod_mean(kk) = Delta_xData_mod_mean(kk-1) + (LiveData_yBar_Diff2(kk+m(j)-1)-LiveData_yBar_Diff2(kk-1))/m(j);
            end
            MHDEV_Output(j,1) = sqrt((mean(Delta_xData_mod_mean.^2)/((Out_tau(j))^2))/6);

            % Calculate the error bar.    
            switch Noise_type
                case "Auto"
                    if floor(N/m(j)) >= 32            
                        log_slope_alpha_approximate = Lag1_NoiseID(LiveData, 1, 2, Phase_or_Fre, m(j)); 
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
            edf_mhdev = EDF_Numerical_sim(log_slope_alpha_approximate,3,m(j), 1, m(j), N);
            Out_ErrBar(j, 1:2) = MHDEV_Output(j,1).*[sqrt(edf_mhdev/chi2inv(0.5+Conf_interval/2, edf_mhdev)), sqrt(edf_mhdev/chi2inv(0.5-Conf_interval/2, edf_mhdev))];
        end
        
    else
        error("Matlab::The data type name is not permitted, 'Frequency' or 'Phase' instead.")
    end
end

function [alpha, d] = Lag1_NoiseID(DataVector, dmin, dmax, DataType, m)

% Input:
    % DataVector :  Imput vector of phase or frequency data.
    % dmin : the minimum order of differencing(default = 0).
    % dmax : the maximum order of differencing.
    % m : Overlapping factor.
% Output:
    % alpha: an estimate of the α of the dominant power law noise type
    % d: (optionally) the value of d.
% Lag 1 autocorrelation. Edit by Tangyi Liao
% Update:2024-11-26,16:11 release: 2.0
% Reference: W. J. Riley and C. A. Greenhall, "Power law noise identification using the lag 1 autocorrelation," 2004 18th (EFTF 2004), Guildford, 2004, pp. 576-580.
    
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

function edf_out = EDF_Numerical_sim(alpha, d, m, F, S, N)
% alpha:frequecy noise exponent
% d： d = 1: flrst-difierence variance (included for completeness)
    % d = 2: Allan variance
    % d = 3: Hadamard variance
% m: m=tau/tau_0;
% F:filter factor
    % F=1:modified variance
    % F=m:unmodified variance
% S:Stride factor
    % S = 1: nonoverlapped estimator
    % S = m: overlapped estimator
% N: number of phase data with sample period tau_0
% edf_out: edf

% Reference: 
% [1] GREENHALL C A, RILEY W J. Uncertainty of Stability Variances Based on Finite Differences, F, 2004 [C].

% Edit by:Tangyin Liao
% First edit : 2022-7-25 21:02 V1.0
% Update : 2024-11-27 11:17 V2.0 , 
    % Fixed some bugs.
    % Add full version output.
    if nargin < 6 
        error('Matlab::Parameter not enough for this function!')
    end
    if alpha+2*d < 1
        error('Matlab::Does not satisfy the restrictions. Should be: alpha + 2d>1')
    end
    
    L = m/F +m*d;
    M = 1+floor(S*(N-L)/m);
    J = min(M,(d+1)*S);
    Jmax = 100; % only in the full version. Suggest value by the reference artical.
    if N < L
        error('Matlab::There are not enough data. Please reduce the value of m.')
    end
    % Main procedure, simplified version
    % edf = (s_z(0,F,alpha,d)^2)*M/BasicSum(J,M,S,F,alpha,d); 

    % Main procedure, full version
    r = M/S;
    table1 = [2/3, 1/3, 7/9, 1/2, 22/25, 2/3; ...
              0.840, 0.345, 0.997, 0.616, 1.141, 0.843;...
              1.079, 0.368, 1.033, 0.607, 1.184, 0.848;...
              NaN, NaN, 1.048, 0.534, 1.180, 0.816;...
              NaN, NaN, 1.302, 0.535, 1.175, 0.777;...
              NaN, NaN, NaN,   NaN,   1.194, 0.703;...
              NaN, NaN, NaN,   NaN,   1.489, 0.702];
    table2 = [1.5, 0.5, 35/18, 1, 2.31, 1.5; ...
              78.6, 25.2, 790, 410, 9950, 6520;...
              2/3, 1/6, 2/3, 1/3, 7/9, 0.5;...
              NaN, NaN, 0.852, 0.375, 0.997, 0.617;...
              NaN, NaN, 1.079, 0.368, 1.033, 0.607;...
              NaN, NaN, NaN,   NaN,   1.053, 0.553;...
              NaN, NaN, NaN,   NaN,   1.302, 0.535];
    table3 = [6, 4, 15.23, 12, 47.8, 40];
    if F == 1 % Modified variances: F = 1, all alpha. This case also applies to unmodified variances when F = m = 1.
        if J<Jmax
            edf_out = (s_z(0,1,alpha,d)^2)*M/BasicSum(J,M,S,1,alpha,d);
        elseif J>Jmax && r>=d+1
            a0 = table1(3-alpha, 2*d-1);
            a1 = table1(3-alpha, 2*d);
            edf_out = r/(a0-a1/r);
        else 
            m_new = Jmax/r;
            edf_out = (s_z(0,1,alpha,d)^2)*Jmax/BasicSum(Jmax,Jmax,m_new,1,alpha,d);
        end
    elseif F>1 && alpha<=0 % Unmodified variances, WHFM to RRFM: F = m, alpha <= 0
        if J<=Jmax
            if m*(d+1)<=Jmax
                edf_out = (s_z(0,m,alpha,d)^2)*M/BasicSum(J,M,S,m,alpha,d);
            else
                edf_out = (s_z(0,1e20,alpha,d)^2)*M/BasicSum(J,M,S,1e20,alpha,d);
            end
        elseif J>Jmax && r>=d+1
            a0 = table2(3-alpha, 2*d-1);
            a1 = table2(3-alpha, 2*d);
            edf_out = r/(a0-a1/r);
        else
            m_new = Jmax/r;
            edf_out = (s_z(0,1e20,alpha,d)^2)*Jmax/BasicSum(Jmax,Jmax,m_new,1e20,alpha,d);
        end
    elseif F>1 && alpha == 1 % Unmodified variances, FLPM: F = m, alpha = 1
        if J<=Jmax         
            edf_out = (s_z(0,m,alpha,d)^2)*M/BasicSum(J,M,S,m,alpha,d);
            % Remark: For this case, m must be less than about 106 to avoid roundo® error.
        elseif J>Jmax && r>=d+1
            a0 = table2(3-alpha, 2*d-1);
            a1 = table2(3-alpha, 2*d);
            b0 = table3(1,2*d-1);
            b1 = table3(1,2*d);
            edf_out = ((b0+b1*log(m))^2)*r/(a0-a1/r);
        else
            m_new = Jmax/r;
            b0 = table3(1,2*d-1);
            b1 = table3(1,2*d);
            edf_out = ((b0+b1*log(m))^2)*Jmax/BasicSum(Jmax,Jmax,m_new,m_new,alpha,d);
        end
    elseif F>1 && alpha == 2 % Unmodified variances, WHPM: F = m, alpha = 2
        K = ceil(r);
        if K<=d
            aa = [];
            for k=1:1:K-1
                aa = aa + (1-k/r)*(factorial(2*d)/(factorial(k)*factorial(d-k)));
            end
            edf_out = M/(1+2*((factorial(d)*factorial(d)/factorial(2*d))^2)*aa );
        else
            a0 = table2(3-alpha, 2*d-1);
            a1 = table2(3-alpha, 2*d);
            edf_out = M/(a0-a1/r);
        end
    end
end

function out_s_w = s_w(t,alpha)
    if t ~= 0
        if alpha == 2
            out_s_w = -abs(t);
        elseif alpha == 1
            out_s_w = (t^2)*log(abs(t));
        elseif alpha == 0
            out_s_w = abs(t^3);
        elseif alpha == -1
            out_s_w = -(t^4)*log(abs(t));
        elseif alpha == -2
            out_s_w = -(abs(t))^5 ;
        elseif alpha == -3
            out_s_w = t^6*log(abs(t));
        elseif alpha == -4
            out_s_w = abs(t)^7;
        end
    else 
        out_s_w = 0;
    end
end

function out_s_x = s_x(t,F,alpha)
    if F < 1e20
        out_s_x = F^2*(2*s_w(t,alpha)-s_w(t-1/F,alpha)-s_w(t+1/F,alpha));
    else
        out_s_x = s_w(t,alpha+2);
    end
end


function out_s_z = s_z(t,F,alpha,d)
    if d == 1
        out_s_z = 2*s_x(t,F,alpha)-s_x(t-1,F,alpha)-s_x(t+1,F,alpha);
    elseif d==2
        out_s_z = 6*s_x(t,F,alpha)-4*s_x(t-1,F,alpha)-4*s_x(t+1,F,alpha)+s_x(t-2,F,alpha)+s_x(t+2,F,alpha);
    elseif d==3
        out_s_z = 20*s_x(t,F,alpha)-15*s_x(t-1,F,alpha)-15*s_x(t+1,F,alpha)+6*s_x(t-2,F,alpha)+6*s_x(t+2,F,alpha)-s_x(t-3,F,alpha)-s_x(t+3,F,alpha);
    end

end

function basicsum = BasicSum(J,M,S,F,alpha,d)
    aa=0;
    for j = 1:1:J-1
        aa = aa + 2*(1-j/M)*(s_z(j/S,F,alpha,d)^2);
    end
    basicsum = s_z(0,F,alpha,d)^2 + (1-J/M)*(s_z(J/S,F,alpha,d)^2) + aa;
end


