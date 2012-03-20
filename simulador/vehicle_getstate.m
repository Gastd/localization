function [vehiclestate] = vehicle_getstate(t,trajectory_name)

vehiclestate = struct('x',0,'y',0,'z',0,'pitch',0,'yaw',0,'roll',0,'dx_dt',0,'dy_dt',0,'dz_dt',0,'d2x_dt2',0,'d2y_dt2',0,'d2z_dt2',0,'dpitch_dt',0,'dyaw_dt',0,'droll_dt',0);
flagtrajectorydefined = 0;

if strcmp(trajectory_name,'trajectory_helix')
   %%% Helix trajectory
   flagtrajectorydefined = 1;
   
   vz = 0.5; % vertical speed in m/s (up is positive)
   wxy = 1*pi/180; % rotational speed in rad/s
   r = 400; % circle radius
   
   vehiclestate.x = r*cos(wxy*t); % North
   vehiclestate.y = r*sin(wxy*t); % East
   vehiclestate.z = -vz*t; % Down
   vehiclestate.dx_dt = -r*wxy*sin(wxy*t);
   vehiclestate.dy_dt =  r*wxy*cos(wxy*t);
   vehiclestate.dz_dt = -vz;
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

if strcmp(trajectory_name,'trajectory_figure8')
    flagtrajectorydefined = 1;
    TrackLength = 2000.0;
    Speed = 25.0;
    CrossOverHeight = 100.0;
    %
    % Grewal & Andrews, Kalman Filtering: Theory and Practice Using MATLAB,
    % John Wiley & Sons, 2008
    %
    %
    % Simulator for vehicle on a figure-eight track, in locally level
    % coordinates, with near-critical banking
    %
    % INPUTS
    %   Time            time in seconds from crossover point
    %   TrackLength     track length in meters
    %   Speed           average vehicle speed in meters per second
    %   CrossOverHeight height of bridge at crossover
    %
    % OUTPUTS
    %   VehState        vector of vehicle state variables:
    %       PosN        Northing from crossover [meters]
    %       PosE        Easting from crossover [meters]
    %       PosD        Downward position wrt median crossover [meters]
    %       VelN        North velocity [m/s]
    %       VelE        East velocity [m/s]
    %       VelD        Downward velocity [m/s]
    %       AccN        North acceleration [m/s/s]
    %       AccE        East acceleration [m/s/s]
    %       AccD        Downward acceleration [m/s/s] (not including gravity)
    %       Roll        Vehicle roll angle [rad]
    %       Pitch       Vehicle pitch angle, up from horizontal [rad]
    %       Heading     Vehicle heading measured clockwise from north [rad]
    %       RollRate    Vehicle rotation rate about its roll axis [rad/s]
    %       PitchRate   Vehicle rotation rate about its pitch axis [rad/s]
    %       YawRate     Vehicle rotation rate about its yaw axis [rad/s]
    %       AccR        Acceleration along vehicle roll axis [m/s/s]
    %       AccP        Acceleration along vehicle pitch axis [m/s/s]
    %       AccY        Acceleration along vehicle yaw axis [m/s/s]
    %                   (not including gravity)
    %       AccR        Vehicle acceleration along its roll axis [m/s/s];
    %       AccP        Vehicle acceleration along its pitch axis [m/s/s];
    %       AccY        Vehicle acceleration along its yaw axis [m/s/s];
    %
    S         = TrackLength/14.94375529901562;  % track scaling parameter
    h         = -CrossOverHeight;
    omega     = 2*pi*Speed/TrackLength;
    G         = 9.8;
    MaxRoll   = asin(3*omega^2*S/G);
    theta     = omega*t;
    s1        = sin(theta);
    c1        = cos(theta);
    s2        = sin(2*theta);
    c2        = cos(2*theta);
    % Northing from crossover [meters]
    PosN      = 3*S*s1;
    % Easting from crossover [meters]
    PosE      = S*s2;
    % Downward position wrt median crossover [meters]
    PosD      = -h*c1/2;
    % North velocity [m/s]
    VelN      = 3*S*c1*omega;
    % East velocity [m/s]
    VelE      = 2*S*c2*omega;
    % Downward velocity [m/s]
    VelD      = h*s1*omega/2;
    % North acceleration [m/s/s]
    AccN      = -3*S*s1*omega^2;
    % East acceleration [m/s/s]
    AccE      = -4*S*s2*omega^2;
    % Downward acceleration [m/s/s] (not including gravity)
    AccD      = h*c1*omega^2/2;
    % Vehicle roll angle [rad]
    Roll      = MaxRoll*s1;
    % Vehicle pitch angle, up from horizontal [rad]
    Pitch     = atan2(h*s1/2,S*(9*c1^2 + 4*c2^2)^(1/2));
    % Vehicle heading measured clockwise from north [rad]
    Heading   = atan2(2*c2,3*c1);

    RateRoll  = MaxRoll*c1*omega;
    % Rate of increase in pitch angle [rad/s]
    RatePitch = 2*h*omega*S*(4*c1*c2^2 + 9*c1 + 8*s1*c2*s2)/(9*c1^2 + 4*c2^2)^(1/2)/(36*S^2*c1^2 + 16*S^2*c2^2 + h^2*(1 - c1^2));
    % Rate of change in vehicle heading [rad/s]
    RateAz    = 6*omega*(s1*c2 - 2*c1*s2)/(9*c1^2 + 4*c2^2);
    
    vehiclestate.x = PosN;
    vehiclestate.y = PosE;
    vehiclestate.z = PosD;
    vehiclestate.dx_dt = VelN;
    vehiclestate.dy_dt = VelE;
    vehiclestate.dz_dt = VelD;
    vehiclestate.d2x_dt2 = AccN;
    vehiclestate.d2y_dt2 = AccE;
    vehiclestate.d2z_dt2 = AccD;
    vehiclestate.pitch = Pitch;
    vehiclestate.yaw = Heading;
    vehiclestate.roll = Roll;
    vehiclestate.dpitch_dt = RatePitch;
    vehiclestate.dyaw_dt = RateAz;
    vehiclestate.droll_dt = RateRoll;
end

if flagtrajectorydefined == 0,
   error(sprintf('Trajectory name %s not defined',trajectory_name));
end
