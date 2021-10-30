%--------------------------------------------------------------------------
%   Single Point Positioning using  pseudoranges & carrier phase measurements 
%  ------------------------------------------------------------------------ 
%  Coder : Doha HASSAN
%  Date  : 26 - 10 - 2021 
%==========================================================================
clear all;
%--------------------------------------------------------------------------
% initialized vectors and matrices
x_satelite = zeros(8,6); dt_satelite = zeros(8,1); Ranges = zeros(8,1); elevation_angle = Ranges; azimuth = Ranges; 
%--------------------------------------------------------------------------
% constants
light_speed = 299792458;         
Wedot= 7.2921151467e-5;     % WGS 84 value of earth’s rotation rate (r/s).
mu = 3.986005e+14;           % WGS 84 value of earth's universal gravitation constant (m^3/s^2)
F = -4.442807633e-10;        % Relativistic correction term constant.
%--------------------------------------------------------------------------
% Input Data
navigation_data  ='eph.dat';
observations  ='rcvr.dat';
unknowns = 4; 
res_thr = 10;
cutoff_thr = 10;
%--------------------------------------------------------------------------
% Output Results
satelite_position_file = 'satelite_position.txt'; file_satelite_position = fopen(satelite_position_file,'w');
receiver_position_file = 'receiver_position.txt'; file_receiver_position = fopen(receiver_position_file,'w');
satelite_residues_file = 'satelite_residues.txt'; file_satelite_residues = fopen(satelite_residues_file,'w');
fprintf(file_satelite_position,'PRN              X(m)                     Y(m)                    Z(m)                    Sat_clock(sec)  \n');
fprintf(file_receiver_position,'Iter.     Delta_X(m)        Delta_Y(m)          Delta_Z(m)      Delta_clock(m)     Final_X(m)      Final_Y(m)     Final_Z(m)    Final_clock(m)');
fprintf(file_receiver_position,'   Sig_X(m)  Sig_Y(m)  Sig_Z(m)  Sig_clock(m) DOP\n');
fprintf(file_satelite_residues,'Iter.      Sat.     Azimuth      elevation_angle     Code RES.(m) \n');
 
%--------------------------------------------------------------------------
% satellite navigation data
[TOW_orb,Sat_PRN_ORB,TOE,a0,a1,a2,e,roota,dn,m_0,dot,i_0,omega,omega_0,omega_dot,cus,cuc,cis,cic,crs,crc] = Broadcast_information(navigation_data);
%--------------------------------------------------------------------------
% observation data
[Time_Of_Week,SV_id, Pseudocode, carrier_phase] = read_measurements(observations);
n = length(Pseudocode);
%--------------------------------------------------------------------------
% compute satellite positions and clocks 
for i = 1:length(SV_id)
    PRN_i = SV_id(i);
    I = find(Sat_PRN_ORB == PRN_i);
    [x_satelite(i,1:3), dt_satelite(i,1)] = Sat_Position((Time_Of_Week(I,1)),Pseudocode(I,1),TOE(I,1),...
        a0(I,1),a1(I,1),a2(I,1), TOE(I,1), m_0(I,1),roota(I,1)^2,...
        dn(I,1),e(I,1),omega(I,1), omega_0(I,1), omega_dot(I,1),i_0(I,1),...
        dot(I,1),cus(I,1),cuc(I,1),cis(I,1),cic(I,1),crs(I,1),crc(I,1),TOE(I,1));
end
%-------------------------------------------------------------------------
% print the satellite positions & clocks 
for i=1:length(SV_id)
    PRN_i = SV_id(i);
    fprintf(file_satelite_position,'G%02d       %+15.3f          %+15.3f          %+15.3f          %+20.15f\n',PRN_i, x_satelite(i,1),x_satelite(i,2),x_satelite(i,3),dt_satelite(i,1));   
end
fclose(file_satelite_position);
%--------------------------------------------------------------------------
% initial receiver coordinates/receiver clock offset
receiver_position   = [-2694685.473; -4293642.366; 3857878.924];
clk_bias = 0;
%--------------------------------------------------------------------------
% geodetic coordinates
[w_latitude, w_longitude, w_altitude] = Wgsxyz2lla(receiver_position);
%--------------------------------------------------------------------------
% loop over code observations
iteration = 1;
while (1)
    while(1)
        count = 1; A = []; calculated_range = []; observed_range = []; P = [];
        azimuth = [];  elevation_angle=[]; Ranges=[];
        for Satelite_index = 1:length(SV_id)
%--------------------------------------------------------------------------
% elevation_angle/azimuth/distance computation
            [azimuth(Satelite_index,1),elevation_angle(Satelite_index,1),Ranges(Satelite_index,1)] = SV_Elevation_Azimuth(receiver_position(1),receiver_position(2),receiver_position(3),w_latitude*pi/180,w_longitude*pi/180,x_satelite(Satelite_index,1),x_satelite(Satelite_index,2),x_satelite(Satelite_index,3));
%------------------------------------------------------------------------
% design of coefficient matrix
            if (elevation_angle(Satelite_index,1) > cutoff_thr)
                A(count,1) = (receiver_position(1) - x_satelite(Satelite_index,1))/Ranges(Satelite_index,1); 
                A(count,2) = (receiver_position(2) - x_satelite(Satelite_index,2))/Ranges(Satelite_index,1); 
                A(count,3) = (receiver_position(3) - x_satelite(Satelite_index,3))/Ranges(Satelite_index,1); 
                A(count,4) =  1;                                       
 
% -------------------------------------------------------------------------
% calculated_range vector
                calculated_range(count,1) = Ranges(Satelite_index,1) + clk_bias - light_speed*dt_satelite(Satelite_index,1); 
% ------------------------------------------------------------------------
% observed range vector
                observed_range(count,1) = Pseudocode(Satelite_index,1) ;
% ------------------------------------------------------------------------
% observation covariance matrix
                P(count, count) = 1 ;    
                count = count+1;     
            end
        end        
% ------------------------------------------------------------------------
% least squares adjustment 
        W = diag((diag(P).^-1));
% normal matrix
        N = (A'*(W)*A);
        if (length(SV_id)>= unknowns) 
% Estimation of unkonwn parameters from Least squares solution
            X   = (N^-1)*A'*(W)*(observed_range - calculated_range);
% Residuals
            V = (A*X + calculated_range)-observed_range;        
        end
        
% Which resilduals are rejected!
        [V_max, reject_V, Sat_PRN_Rej] = Max_Residue(V, SV_id);
        if (V_max > res_thr)
            [SV_id,x_satelite,dt_satelite,Pseudocode] = Reject_code_measurements(SV_id,x_satelite,dt_satelite,Pseudocode,Sat_PRN_Rej);
        else
            break;
        end
    end
    
% Solution update
    receiver_position   = receiver_position + X(1:3);
    clk_bias = clk_bias + X(4);
        
% Convergence check
    delta_x = norm(X(1:3));   
        
% a posteriori covariance matrix
    posteriori_cov = (V'*(W)*V)/(count - unknowns);
    C_xx = posteriori_cov*(N^-1);
    cov_XR  = C_xx(1:3,1:3);
    var_dtR = C_xx(4,4);
    std_x  = sqrt(cov_XR(1,1));  std_y  = sqrt(cov_XR(2,2));   std_z  = sqrt(cov_XR(3,3));
    std_dtR = sqrt(var_dtR);
% ------------------------------------------------------------------------        
    % DOP
    cov_XYZ = (A'*A)^-1;
    PDOP   = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
% ------------------------------------------------------------------------        
    % print position results
    fprintf(file_receiver_position,'%d        %+12.5f  %+12.5f  %+12.5f  %+11.3f   ',iteration,X(1),X(2),X(3),X(4));   
    fprintf(file_receiver_position,'%+12.3f    %+12.3f    %+12.3f    %+12.3f',receiver_position(1),receiver_position(2),receiver_position(3),clk_bias);  
    fprintf(file_receiver_position,'    %+8.3f  %+8.3f  %+8.3f  %+8.3f     %+5.2f\n',std_x,std_y,std_z,std_dtR, PDOP);
% ------------------------------------------------------------------------        
    % print residues
    for j=1:length(SV_id)
        Sat_PRN_p = SV_id(j); I=find(SV_id==Sat_PRN_p);
        if (azimuth(I)< 0), azimuth(I)=azimuth(I)+180; end
        fprintf(file_satelite_residues,'%d              G%02d     %6.2f      %5.2f         %+6.2f  \n',iteration,Sat_PRN_p,azimuth(I),elevation_angle(I),V(I));
    end
    fprintf(file_satelite_residues,'                                           \n');
% ------------------------------------------------------------------------        
    % Convergence check
    if (delta_x < 1e-4) 
        break; 
    end
    iteration = iteration+1;
    if (iteration > 10)
        break; 
    end  
end
fclose(file_receiver_position); fclose(file_satelite_residues);
fclose('all');
% ------------------------------------------------------------------------        

