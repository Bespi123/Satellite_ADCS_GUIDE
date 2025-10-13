function parameters = simulationParameters(~)
    % --- CubeSat Physical Parameters ---
    % Inertia tensor of the CubeSat [kg*m^2].
    parameters.initialValues.I = 1E-3 * [3.084 0.082 -0.054;
                                         0.082 3.132 0.016;
                                        -0.054 0.016 3.540];
    % Initial attitude quaternion [w, x, y, z]'.
    parameters.initialValues.q0 = [1; 0; 0; 0];
    % Initial angular velocity [rad/s].
    parameters.initialValues.Wo = [0; 0; 0];
    
    % --- System Timing and Sampling Parameters ---
    % Sample time for the gyroscope sensor [s].
    parameters.sensors.gyro.Ts      = 0.005; % 5 ms
    % Sample time for other attitude sensors (Accel, Mag, Star Tracker) [s].
    parameters.sensors.attitude.Ts  = 0.01;  % 10 ms
    % Sample time for the control loop execution [s].
    parameters.control.Ts           = 0.005; % 5 ms
    
    % --- Sensor Model Parameters ---

    %%% Accelerometer sensor
    % Gravitational acceleration vector in inertial frame 
    parameters.sensors.acc.g_I  = [0,0,1]';  
    % Accelerometer constant bias [m/s^2].
    parameters.sensors.acc.bias = [0.013127156; 0.010186677; 0.000138897];
    % Accelerometer noise standard deviation [m/s^2].
    parameters.sensors.acc.std  = [0.011817097; 0.012306208; 0.015949601];
    
    % Magnetometer reference direction (in inertial frame)
    parameters.sensors.mag.m_I  = [0.7071 ,0 ,0.7071]';  
    % Magnetometer constant bias [uT or normalized units].
    parameters.sensors.mag.bias = [0; 0; 0];
    % Magnetometer noise standard deviation [uT or normalized units].
    parameters.sensors.mag.std  = [0.011817097; 0.012306208; 0.015949601];
    
    % Gyroscope constant bias [rad/s].
    parameters.sensors.gyro.bias = [-0.0000541; -0.0001553; -0.00006458];
    % Gyroscope noise standard deviation (Angular Random Walk) [rad/s].
    parameters.sensors.gyro.std  = [0.009; 0.008; 0.007];
    %%% Gyro filter time constant
    parameters.sensors.gyro.tau = 0.1;
    
    % Star Tracker constant bias.
    parameters.sensors.star.bias = 1E-6 * [1; 1; 1];
    % Star Tracker noise standard deviation.
    parameters.sensors.star.std  = 1E-3 * [1; 1; 1];
    % Flag to enable (1) or disable (0) the star tracker.
    parameters.sensors.star.enable = 1;
    % Number of stars to be simulated by the star tracker.
    parameters.sensors.star.numberOfStars = 10;

    % --- AHRS Selection ---
    % Selector: 1 for MADGWICK, 2 for EKF, 3 for UKF.
    parameters.ahrs.flag = 'EKF';
    % Flag: 1 to enable AHRS feedback to controller, 0 to use ideal state.
    parameters.ahrs.enable = 1;

    % MADGWICK parameters
    parameters.ahrs.madwick.beta = 1;
    parameters.ahrs.madwick.bias_gain = 0.01;

    % EKF parameters
    % Flag: 1 to use sensor model noise for EKF, 0 to use custom values.
    parameters.ahrs.ekf.equalModel = 0;
    % EKF accelerometer noise standard deviation.
    parameters.ahrs.ekf.acc_std    = [0.011817097; 0.012306208; 0.015949601];
    % EKF magnetometer noise standard deviation.
    parameters.ahrs.ekf.mag_std    = [0.011817097; 0.012306208; 0.015949601];
    % EKF gyroscope noise standard deviation.
    parameters.ahrs.ekf.gyro_std   = [0.009; 0.008; 0.007];
    % EKF star tracker noise standard deviation.
    parameters.ahrs.ekf.star_std   = 1E-3 * [1; 1; 1];
    
    % UKF parameters
    parameters.ahrs.ukf.alpha = 1e-3;
    parameters.ahrs.ukf.beta  = 2.0;
    parameters.ahrs.ukf.kappa = 0.0;
    parameters.ahrs.ukf.gyro_std = [0.009; 0.008; 0.007];
    parameters.ahrs.ukf.acc_std  = [0.011817097; 0.012306208; 0.015949601];
    parameters.ahrs.ukf.mag_std  = [0.011817097; 0.012306208; 0.015949601];
    parameters.ahrs.ukf.star_std  = 1E-3 * [1; 1; 1];

    % --- Controller Setpoint Parameters ---
    % Desired target attitude quaternion [w, x, y, z]'.
    parameters.setPoint.qd = [0.9808; 0.1691; 0.0933; 0.0277];
    % Desired target angular velocity [rad/s].
    parameters.setPoint.Wd = [0; 0; 0];

    % --- Controller Selection ---
    % Selector: 1 for Feedback Controller, 2 for Boskovic's Controller.
    parameters.controller.selector = 1;

    % --- Feedback Controller Parameters ---
    % Derivative gain matrix.
    parameters.controller.feedback.Keye = 0.2 * eye(3);
    % Proportional gain matrix.
    parameters.controller.feedback.Peye = 0.3 * eye(3);

    % --- Boskovic Controller Parameters ---
    parameters.controller.boskController.delta = 0.5;
    parameters.controller.boskController.gamma = 0.001;
    parameters.controller.boskController.Umax  = 0.10;
    parameters.controller.boskController.k0    = 1;

    % --- Disturbance Torque Parameters ---
    % Constant disturbance torque components [Nm].
    parameters.disturbances.constant.Tdx = 0.001;
    parameters.disturbances.constant.Tdy = 0.001;
    parameters.disturbances.constant.Tdz = 0.001;
    % Flag to enable (2) or disable (1) constant disturbances.
    parameters.disturbances.constant.enable = 1;
    
    % Parameters for variable (sinusoidal) disturbances [Nm].
    parameters.disturbances.variable.A   = 0.01;
    parameters.disturbances.variable.w1  = 1.00;
    parameters.disturbances.variable.w2  = 2.00;
    parameters.disturbances.variable.w3  = 3.00;
    % Flag to enable (2) or disable (1) variable disturbances.
    parameters.disturbances.variable.enable = 1;

    % --- Reaction Wheels Configuration ---
    % Number of reaction wheels in the assembly.
    parameters.rw.number = 3;
    % Maximum momentum storage per wheel [Nms].
    parameters.rw.maxH   = 1;
    % Mounting angles for wheels (example for a 3-wheel orthogonal config).
    parameters.rw.beta1  = 0;
    parameters.rw.beta2  = 0;
    % Distribution matrix for mapping control torque to wheel torques.
    parameters.rw.W      = eye(3);

    % --- Reaction Wheel DC Motor Model ---
    parameters.rw.motor.kt   = 0.0254;   % Torque constant [Nm/A].
    parameters.rw.motor.Jrw  = 3.475E-5; % Rotor inertia [kg*m^2].
    parameters.rw.motor.b    = 5.58E-4;  % Viscous friction coefficient [Nms].
    parameters.rw.motor.c    = 0.0000;   % Coulomb friction [Nm].
    parameters.rw.motor.L    = 0.000403; % Armature inductance [H].
    parameters.rw.motor.R    = 1.16;     % Armature resistance [Ohms].
    parameters.rw.motor.ke   = 0.0254;   % Back-EMF constant [Vs/rad].
    
    % --- Reaction Wheel DC initial conditions---
    parameters.rw.motor.init.w_rw = 0;
    parameters.rw.motor.init.current = 0; 
    
    % --- Reaction Wheel PID Speed Controller Gains ---
    parameters.rw.motor.pid.kp = 0.504519552012394;
    parameters.rw.motor.pid.ki = 28.7345189644697;
    parameters.rw.motor.pid.kd = -0.000554999924740984;
    parameters.rw.motor.pid.Nu = 289.873462888377; % Derivative filter coefficient.
end