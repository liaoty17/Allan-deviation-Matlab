function edf_out = EDF_Numerical_sim(alpha, d, m, F, S, N)
% Using for estimate the noise evidence factor.

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


