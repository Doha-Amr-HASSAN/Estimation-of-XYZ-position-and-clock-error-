%   This function computes the satellite position and satellite clock coreection from the broadcast ehemeris.
%  Coder : Doha HASSAN
%  Date  : 26 - 10 - 2021
%--------------------------------------------------------------------------
function [Satelite_pos, Satelite_Clc_Corr] = Sat_Position(time,pr,toc,f0,f1,f2, toe, M0,A,deln,e,omg,OMG0,OMGd,i0,idot,cus,cuc,cis,cic,crs,crc,TOEs)

Light_Speed = 299792458.d0;
time = time-(pr/Light_Speed);
delta_t = eph2clk(time, toc, f0, f1, f2);
time = time - delta_t;
[Satelite_pos, Satelite_Clc_Corr] = eph2pos(time,toe,M0,A,deln,e,omg,OMG0,OMGd,i0,idot,cus,cuc,cis,cic,crs,crc,TOEs,f0,f1,f2);
end 
function eph2clk = eph2clk(time, toc, f0, f1, f2)
t=time-toc;
for i=1:2
    t=t-(f0+f1*t+f2*t*t);
end
eph2clk=f0+f1*t+f2*t*t;
end
function [Satelite_pos, Satelite_Clc_Corr] = eph2pos(time,TOE,M0,A,deln,e,omg,OMG0,OMGd,i0,idot,cus,cuc,cis,cic,crs,crc,toes,f0,f1,f2)
Light_Speed = 299792458.d0;
MU_GPS = 3.9860050d14;          % gravitational constant
OMGA_E   = 7.292115d-5;           % earth angular velocity (rad/s) ref (2) 
RTOL_KEPLER = 1d-13;            % relative tolerance for Kepler equation 
Kepler_Max_Iter = 30;           % max number of iteration of Kelpler 

if(A<=0.d0)
    Satelite_pos = 0.d0; Satelite_Clc_Corr = 0.d0;
    return;
end

tk=(time-TOE);
mu=MU_GPS; omge1 = OMGA_E;
M=M0+(sqrt(mu/(A*A*A))+deln)*tk;
n=0;E=M;Ek=0.d0;

while(abs(E-Ek)>RTOL_KEPLER && n < Kepler_Max_Iter)
    Ek=E; E=E-(E-e*sin(E)-M)/(1.d0-e*cos(E));
    n=n+1;
end

if(n >= Kepler_Max_Iter), return; end
sinE=sin(E); cosE=cos(E);

u=atan2(sqrt(1.d0-e*e)*sinE,cosE-e)+omg;
r=A*(1.d0-e*cosE);
i=i0+idot*tk;
sin2u=sin(2.d0*u); cos2u=cos(2.d0*u);
u=u+cus*sin2u+cuc*cos2u;
r=r+crs*sin2u+crc*cos2u;
i=i+cis*sin2u+cic*cos2u;
x=r*cos(u); y=r*sin(u); cosi=cos(i);


O=OMG0+(OMGd-omge1)*tk-omge1*toes;
sinO=sin(O); cosO=cos(O);
Satelite_pos(1)=x*cosO-y*cosi*sinO;
Satelite_pos(2)=x*sinO+y*cosi*cosO;
Satelite_pos(3)=y*sin(i);

tk=(time-TOE);
Satelite_Clc_Corr = f0+f1*tk+f2*tk*tk;

% relativity correction 
Satelite_Clc_Corr = Satelite_Clc_Corr-2.d0*sqrt(mu*A)*e*sinE/(Light_Speed*Light_Speed);
end 
%--------------------------------------------------------------------------
