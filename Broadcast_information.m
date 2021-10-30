% This function reads the information for the satellite positions and satellite clock correction.
%  Coder : Doha HASSAN
%  Date  : 26 - 10 - 2021 
%--------------------------------------------------------------------------
function [TOW,Sat_PRN_ORB,TOE,a0,a1,a2,e,roota,dn,m0,dot,i0,omega,omega0,omega_dot,cus,cuc,cis,cic,crs,crc] = Broadcast_information(eph)
% initialization        
TOW         = zeros(8,1); Sat_PRN_ORB = zeros(8,1);TOE = zeros(8,1);
a0          = zeros(8,1);a1  = zeros(8,1); a2  = zeros(8,1);
e           = zeros(8,1);
roota       = zeros(8,1);
dn          = zeros(8,1);
m0          = zeros(8,1);
dot         = zeros(8,1);
i0          = zeros(8,1);
omega       = zeros(8,1);omega0   =zeros(8,1);omega_dot =zeros(8,1);

cus = zeros(8,1);
cuc = zeros(8,1);
cis = zeros(8,1);
cic = zeros(8,1);
crs = zeros(8,1);
crc = zeros(8,1);

Satelite_index=1;

file_id = fopen(eph,'rt');

while(~feof(file_id))
    
    line=fgetl(file_id);
    if (length(line)>1)
        str=sscanf(line(1:16),'%c');
        TOW(Satelite_index,1)=str2double(str);
        str=sscanf(line(18:19),'%c');
        Sat_PRN_ORB(Satelite_index,1)=str2double(str);
        str=sscanf(line(41:53),'%c');
        TOE(Satelite_index,1)=str2double(str);
        str=sscanf(line(57:70),'%c');
        a0(Satelite_index,1)=str2double(str);
        str=sscanf(line(74:87),'%c');
        a1(Satelite_index,1)=str2double(str);
        str=sscanf(line(92:104),'%c');
        a2(Satelite_index,1)=str2double(str);
        str=sscanf(line(134:146),'%c');
        e(Satelite_index,1)=str2double(str);
        str=sscanf(line(159:171),'%c');
        roota(Satelite_index,1)=str2double(str);
        str=sscanf(line(184:196),'%c');
        dn(Satelite_index,1)=str2double(str);
        str=sscanf(line(208:221),'%c');
        m0(Satelite_index,1)=str2double(str);
        str=sscanf(line(233:246),'%c');
        omega(Satelite_index,1)=str2double(str);
        str=sscanf(line(258:271),'%c');
        omega0(Satelite_index,1)=str2double(str);
        str=sscanf(line(284:296),'%c');
        i0(Satelite_index,1)=str2double(str);
        str=sscanf(line(300:313),'%c');
        omega_dot(Satelite_index,1)=str2double(str);
        str=sscanf(line(317:330),'%c');
        dot(Satelite_index,1)=str2double(str);
        
        str=sscanf(line(335:347),'%c');
        cus(Satelite_index,1)=str2double(str);
        str=sscanf(line(351:364),'%c');
        cuc(Satelite_index,1)=str2double(str);
        str=sscanf(line(368:381),'%c');
        cis(Satelite_index,1)=str2double(str);
        str=sscanf(line(385:398),'%c');
        cic(Satelite_index,1)=str2double(str);
        str=sscanf(line(402:415),'%c');
        crs(Satelite_index,1)=str2double(str);
        str=sscanf(line(420:432),'%c');
        crc(Satelite_index,1)=str2double(str);
        
        Satelite_index=Satelite_index+1;    
    end
end
fclose(file_id);