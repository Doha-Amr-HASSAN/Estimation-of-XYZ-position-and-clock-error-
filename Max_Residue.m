 %  Coder : Doha HASSAN
%  Date  : 26 - 10 - 2021 
%--------------------------------------------------------------------------
function [max_V,reject_V, Sat_PRN_Rej] = Max_Residue(res,iprn_obs)
max_V=0; reject_V=[];

for i=1:length(res)
    if (abs(res(i))>max_V)
       max_V   = abs(res(i));
       reject_V  = i;
       Sat_PRN_Rej = iprn_obs(i);
    end   
end