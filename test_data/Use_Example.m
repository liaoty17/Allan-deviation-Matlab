clc;clear
close all;
addpath("..\")
%% import data
tic
fre_or_phase = 1;
if fre_or_phase == 1
    datain = readmatrix(".\Noise_gene_1127_phas.txt");
    tau0 = 1;
    data_type = "Phase";
else
    datain = readmatrix(".\Noise_gene_1127_freq.txt");
    tau0 = 1;
    data_type = "Frequency";
end
len_m = floor(log2(length(datain)/2)); 
m = 2.^(0:1:len_m-1)'; % Octave: Calculate using the exponential sequence of 2.
tauin = m';

% [Out_tau, DEV_Output, Out_ErrBar, alpha_memory] = AllanTools_Matlab(DeviationType, LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)
[out_tau, dev_out, out_err, alpha] = Allan_ADEV(datain, tau0, data_type, tauin); % Use default values for the remaining parameters.
total_out = [out_err(:,1), dev_out, out_err(:,2)]

toc
figure
errorbar(out_tau, dev_out, dev_out-out_err(:,1), out_err(:,2)-dev_out, "LineWidth", 1.5)
set(gca,"XScale","log","YScale","log","FontSize",12)
xlabel("Average time $\tau$","FontSize",14,"FontWeight","normal","FontName","Times New Roman","Interpreter","latex")
ylabel(strcat("Allan Deviation", {32} ,"$\sigma_{y}(\tau)$"),"Interpreter","latex","FontSize",14,"FontWeight","normal","FontName","Times New Roman")
