function [vehiclestate,lastP,lastY] = vehicle_getstate(t,trajectory_name,lastPitch,lastYaw)

vehiclestate = struct('x',0,'y',0,'z',0,'pitch',0,'yaw',0,'roll',0,'dx_dt',0,'dy_dt',0,'dz_dt',0,'d2x_dt2',0,'d2y_dt2',0,'d2z_dt2',0,'dpitch_dt',0,'dyaw_dt',0,'droll_dt',0);
flagtrajectorydefined = 0;

if strcmp(trajectory_name,'trajectory_helix')
   %%% Helix trajectory
   flagtrajectorydefined = 1;
   
   vz = 0.5; % vertical speed in m/s
   wxy = 1*pi/180; % rotational speed in rad/s
   r = 400; % circle radius
   
   vehiclestate.x = r*cos(wxy*t);
   vehiclestate.y = r*sin(wxy*t);
   vehiclestate.z = vz*t;
   vehiclestate.dx_dt = -r*wxy*sin(wxy*t);
   vehiclestate.dy_dt =  r*wxy*cos(wxy*t);
   vehiclestate.dz_dt = vz;
   vehiclestate.d2x_dt2 = -r*wxy*cos(wxy*t)*wxy;
   vehiclestate.d2y_dt2 = -r*wxy*sin(wxy*t)*wxy;
   vehiclestate.d2z_dt2 = 0;
   
   trajtrangent = [vehiclestate.dx_dt vehiclestate.dy_dt vehiclestate.dz_dt]';
   trajtrangent = trajtrangent / sqrt(trajtrangent'*trajtrangent); % tangent vector to the trajectory (unity length)
   n = cross([1; 0; 0],[0; 1; 0]); % normal to XY plan (unity length)
   vehiclestate.pitch = -asin(n'*trajtrangent); 
   vehiclestate.yaw = atan2(vehiclestate.dy_dt,vehiclestate.dx_dt); 
   vehiclestate.roll = -pi/2*tanh(wxy);
   vehiclestate.dpitch_dt = -1/sqrt(1-(n'*trajtrangent)^2)*n'*[-r*wxy*cos(wxy*t)*wxy; -r*wxy*sin(wxy*t)*wxy; 0];
   vehiclestate.dyaw_dt = wxy;
   vehiclestate.droll_dt = 0;
end

if strcmp(trajectory_name,'trajectory_hovering')
   %%% Hovering trajectory
   flagtrajectorydefined = 1;
  
   vehiclestate.x = 5;
   vehiclestate.y = -5;
   vehiclestate.z = 0;
   vehiclestate.dx_dt = 0;
   vehiclestate.dy_dt = 0;
   vehiclestate.dz_dt = 0;
   vehiclestate.d2x_dt2 = 0;
   vehiclestate.d2y_dt2 = 0;
   vehiclestate.d2z_dt2 = 0;
   
   tr = 100;
   d = rem(t,tr);
   
   angle = (pi/3)*sin(t*2*pi*3/tr);
   dangle_dt = (pi/3)*cos(t*2*pi*3/tr)*2*pi*3/tr;
   
   if d < tr/3
      vehiclestate.pitch = angle;
      vehiclestate.yaw = 0;  
      vehiclestate.roll = 0;
      vehiclestate.dpitch_dt = dangle_dt;
      vehiclestate.dyaw_dt = 0;
      vehiclestate.droll_dt = 0;
   end
   if (d >= tr/3) & (d < 2*tr/3) 
      vehiclestate.pitch = 0;
      vehiclestate.yaw = angle;  
      vehiclestate.roll = 0;
      vehiclestate.dpitch_dt = 0;
      vehiclestate.dyaw_dt = dangle_dt;
      vehiclestate.droll_dt = 0;
   end
   if d >= 2*tr/3
      vehiclestate.pitch = 0;
      vehiclestate.yaw = 0;  
      vehiclestate.roll = angle;
      vehiclestate.dpitch_dt = 0;
      vehiclestate.dyaw_dt = 0;
      vehiclestate.droll_dt = dangle_dt;
   end
end

if flagtrajectorydefined == 0,
   error(sprintf('Trajectory name %s not defined',trajectory_name));
end
