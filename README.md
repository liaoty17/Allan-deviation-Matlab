# Allan-deviation-Matlab
**Allan deviation(Standard, Overlapping, Modified, Time) and Hadamard deviation calculate tools in Matlab.**
If you simply want to analyze your data using this program, it is recommended to use the more portable full version, see the AllanTools_Matlab in the <allan_easy_to_use> folder. If you want to embed the code into your program, you can use the individual functions.

Using example of Allan_ADEV:
[Out_tau, OADEV_Output, Out_ErrBar, alpha_memory] = Allan_OADEV(LiveData, tau0, Phase_or_Fre, tau,  Conf_interval, Noise_type)
Input:
    LiveData : the data using for calculate the allan deviation, should be 1-D.
    tau0 : time interval.
    tau : Average time sequence.
    Phase_or_Fre : input data type. Parameter type, "String"
        "Phase" : Phase data.
        "Frequency" : Frequency data.
    Conf_interval : Confidence interval, using for calculate the error bar.
    Noise_type : the noise type of your data, if in doubt, "Auto" is recommended. Parameter type, "String";
        "RWFM" : Random walk frequency modulation, h_{-2}
        "FFM" : Flicker frequency modulation, h_{-1}
        "WFM" : White frequency modulation, h_{0}
        "FPM" : Flicker phase modulation, h_{1}
        "WPM" : White phase modulation, h_{2}
Output:
    Out_tau: Average time sequence.
    ADEV_Output: Allan deviation(No Overlapped)
    Out_ErrBar: 
        Column 1: Min sigma
        Column 2: max sigma
Reference: Riley, W. and Howe, D. (2008), Handbook of Frequency Stability Analysis, Special Publication (NIST SP)...
