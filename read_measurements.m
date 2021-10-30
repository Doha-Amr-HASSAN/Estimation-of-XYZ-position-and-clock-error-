%   This function reads the code observations and carrier phase measurements 
%  Coder : Doha HASSAN
%  Date  : 26 - 10 - 2021 
%--------------------------------------------------------------------------
function [Time_Of_Week , SV_id,Pseudocode,carrier_phase] = read_measurements(OBSERVATIONS)    
Time_Of_Week = zeros(8,1);
SV_id = zeros(8,1);
Pseudocode = zeros(8,1);
carrier_phase  = zeros(8,1);
fraction   = zeros(8,1);
Satelite_index = 1;
file_id = fopen(OBSERVATIONS,'rt');
while(~feof(file_id))
    LINE = fgetl(file_id);
    if (length(LINE)>1)
        STR                             = sscanf(LINE(1:16),'%c');
        Time_Of_Week(Satelite_index,1)           = str2double(STR);
        STR = sscanf(LINE(18:19),'%c');
        SV_id(Satelite_index,1)         = str2double(STR);
        STR = sscanf(LINE(25:36),'%c');
        Pseudocode(Satelite_index,1)    = str2double(STR);
        STR                             = sscanf(LINE(43:49),'%c');
        carrier_phase(Satelite_index,1) = str2double(STR);
        STR                             = sscanf(LINE(51:54),'%c');
        fraction(Satelite_index,1)      = str2double(STR);
        carrier_phase(Satelite_index,1) = carrier_phase(Satelite_index,1)+(fraction(Satelite_index,1)/2048);
        
        Satelite_index = Satelite_index+1;    
    end
end
fclose(file_id);