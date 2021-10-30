%  Coder : Doha HASSAN
%  Date  : 13-10-2021 
%--------------------------------------------------------------------------
function [azimuth,elevation_angle,geometric_distance] = SV_Elevation_Azimuth(x_reciever,y_reciever,z_reciever,phi,plam,x_satelite,y_satelite,z_satelite)
rec_lat = phi - pi/2.0;
rec_lon = plam - pi;
sr_latitude = sin(rec_lat);
cr_latitude = cos(rec_lat);
sr_longitude = sin(rec_lon); 
cr_longitude = cos(rec_lon);

delta_x = x_satelite-x_reciever;
delta_y = y_satelite-y_reciever;
delta_z = z_satelite-z_reciever;

delta_u = cr_latitude*cr_longitude*delta_x+ cr_latitude*sr_longitude*delta_y -sr_latitude*delta_z;
delta_v = sr_longitude*delta_x      - cr_longitude*delta_y;
delta_w = sr_latitude*cr_longitude*delta_x+ sr_latitude*sr_longitude*delta_y+ cr_latitude*delta_z;

geometric_distance = sqrt(delta_u^2+delta_v^2+delta_w^2);
azimuth  = atan2(delta_v,delta_u)/pi*180;
elevation_angle  = asin(delta_w/geometric_distance)/pi*180;