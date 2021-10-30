% This function removes the rejected code residues
%  Coder : Doha HASSAN
%  Date  : 26 - 10 - 2021 
%--------------------------------------------------------------------------
function [Sat_PRN_ORB,x_satelite,dt_satelite,Pseudocode] = Reject_code_measurements(Sat_PRN_ORB,x_satelite,dt_satelite,Pseudocode,Sat_PRN_ORB_reject)
index = find(Sat_PRN_ORB == Sat_PRN_ORB_reject);
Sat_PRN_ORB(index)=[];
x_satelite(index,:)=[];
dt_satelite(index)=[];
Pseudocode(index)=[];