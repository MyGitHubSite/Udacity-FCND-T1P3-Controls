
 
# Position control gains
kpPosXY = 40
kpPosZ = 55
KiPosZ = 50

# Velocity control gains
kpVelXY = 15
kpVelZ = 15

# Angle control gains
kpBank = 12    
kpYaw = 4     

# Angle rate gains
kpPQR = 40, 40, 10  


VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  
  
float new_mass = mass * 1.25;
  // cmd.desiredThrustsN[0] = new_mass * 9.81f / 4.f;                // * 1.2114; // front left
  // cmd.desiredThrustsN[1] = new_mass * 9.81f / 4.f;                // * 1.2092; // front right
  // cmd.desiredThrustsN[2] = new_mass * 9.81f / 4.f;                // * 1.2140; // rear left
  // cmd.desiredThrustsN[3] = new_mass * 9.81f / 4.f;                // * 1.2140; // rear right

  float l = L / M_SQRT2;					   // perpendicular length to the axes
  
  float F = collThrustCmd;
  float Fx = momentCmd.x / l;
  float Fy = momentCmd.y / l;
  float Fz = momentCmd.z / kappa;

  float F1 = (F + Fx + Fy - Fz) / 4;
  float F2 = (F - Fx + Fy + Fz) / 4;
  float F3 = (F + Fx - Fy + Fz) / 4;
  float F4 = (F - Fx - Fy - Fz) / 4;
  float total_thrust = F1 + F2 + F3 + F4;

  assert(total_thrust = collThrustCmd);
    
  cmd.desiredThrustsN[0] = CONSTRAIN(F1, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[1] = CONSTRAIN(F2, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[2] = CONSTRAIN(F3, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[3] = CONSTRAIN(F4, minMotorThrust, maxMotorThrust);
  
  
  V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  
    V3F Inertia = V3F(Ixx, Iyy, Izz);

  // output to Quadcopter
  momentCmd = Inertia * kpPQR * (pqrCmd - pqr);
  
  
  V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
 
 
  
    // get inputs to roll-pitch controller

  // from altitude controller  
  float collAccelCmd = -collThrustCmd / mass;     // convert thrust u1 (N) to m/s^2

  // from lateral controller
  float x_c_dot_dot = accelCmd.x;
  float y_c_dot_dot = accelCmd.y;

  // from quad
  float b_x = R(0,2);       // roll angle phi
  float b_y = R(1,2);       // roll angle theta
  // float b_z = R(2, 2);   // roll angle psi   note calculated in yaw controller

  //set target p, q angles   note r is calculated by yaw controller
  float b_x_target = atan(x_c_dot_dot / collAccelCmd);   // p angle
  float b_y_target = atan(y_c_dot_dot / collAccelCmd);   // q angle
  
  b_x_target = CONSTRAIN(b_x_target, -maxTiltAngle, maxTiltAngle);
  b_y_target = CONSTRAIN(b_y_target, -maxTiltAngle, maxTiltAngle);

  float b_c_x_dot = kpBank * (b_x_target - b_x);   // pitch rate
  float b_c_y_dot = kpBank * (b_y_target - b_y);   // roll rate

  // convert coordinate frames
  pqrCmd.x = 1.0 / R(2,2) * (R(1,0) * b_c_x_dot - R(0,0) * b_c_y_dot);
  pqrCmd.y = 1.0 / R(2,2) * (R(1,1) * b_c_x_dot - R(0,1) * b_c_y_dot);
  
  
  float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
 
 
    float g = 9.81f;																	// gravitational force m/s^2
  float z = posZ;															// z position
  float z_dot = velZ;														// z velocity
  float z_Cmd = posZCmd;													// z position command
  float z_dot_Cmd = CONSTRAIN(velZCmd, -maxAscentRate, maxDescentRate);     // z velocity command
  float z_dot_dot_Cmd = accelZCmd;											// z acceleration command
  float z_error = z_Cmd - z;												// z position error
  float z_dot_error = z_dot_Cmd - z_dot;									// z velocity error

  integratedAltitudeError += z_error * KiPosZ * dt;

  float P = (z_error * kpPosZ);
  float I = (z_dot_dot_Cmd + integratedAltitudeError);
  float D = (z_dot_error * kpVelZ);

  float z_dot_dot = P + I + D;

  thrust = (g - z_dot_dot) * mass / R(2,2);   // output to Roll-Pitch controller
  
  
  V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  
    V3F pos_error = posCmd - pos;
  V3F vel_error = velCmd - vel;

  accelCmd.x = accelCmd.x + (pos_error.x * kpPosXY) + (vel_error.x * kpVelXY);
  accelCmd.y = accelCmd.y + (pos_error.y * kpPosXY) + (vel_error.y * kpVelXY);

  // outputs to Roll-Pitch controller
  accelCmd.x = CONSTRAIN(accelCmd.x, -maxAccelXY, maxAccelXY);
  accelCmd.y = CONSTRAIN(accelCmd.y, -maxAccelXY, maxAccelXY);
  
  float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  
    float b = 2 * M_PI;
  yaw = fmodf(yaw, b);

  float yaw_error = yawCmd - yaw;

  // output to Body Rate controller
  yawRateCmd = (yaw_error * kpYaw);
  
