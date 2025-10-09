function parameters = simulationParameters(~)
    % This function defines all configurable parameters for the CubeSat ADCS simulation.
    % It returns a single structure 'parameters' that contains all nested sub-structures
    % for different components of the simulation (e.g., initial values, sensors, control).

    % --- CubeSat Physical Parameters ---
    % Inertia tensor of the CubeSat [kg*m^2]. This represents the distribution of mass
    % and how the satellite resists rotational motion. Off-diagonal terms represent
    % products of inertia.
    parameters.initialValues.I = 1E-3 * [3.084 0.082 -0.054;
                                         0.082 3.132 0.016;
                                        -0.054 0.016 3.540];
    % Initial attitude quaternion [w, x, y, z]' (scalar-first format). 
    % Represents the initial orientation of the CubeSat relative to the inertial frame.
    parameters.initialValues.q0 = [1; 0; 0; 0]; % Starts with no rotation.
    
    % Initial angular velocity of the CubeSat body [rad/s].
    parameters.initialValues.Wo = [0; 0; 0]; % Starts stationary.
    
    % --- System Timing and Sampling Parameters ---
    % These parameters simulate a digital system where sensors and control loops
    % run at discrete time intervals, not continuously.
    
    % Sample time for the gyroscope sensor [s]. Gyros are often sampled at a higher
    % rate because they are used for high-frequency state propagation in the AHRS.
    parameters.sensors.gyro.Ts      = 0.005; % 200 Hz
    
    % Sample time for other attitude sensors (Accelerometer, Magnetometer, Star Tracker) [s].
    % These sensors are typically slower and used for correcting the gyro drift.
    parameters.sensors.attitude.Ts  = 0.01;  % 100 Hz
    
    % Sample time for the control loop execution [s]. This defines how often the
    % control law is calculated and a new command is sent to the actuators.
    parameters.control.Ts           = 0.005; % 200 Hz
    
    % --- Sensor Model Parameters ---
    % These sections model the imperfections of real-world sensors.
    
    %%% Accelerometer sensor
    % Gravitational acceleration vector in the inertial frame (e.g., pointing towards Earth).
    % This is the 'true' vector that the accelerometer measures to help determine 'down'.
    parameters.sensors.acc.g_I  = [0,0,1]';  
    % Accelerometer constant bias [m/s^2]. A persistent offset in the measurement.
    parameters.sensors.acc.bias = [0.013127156; 0.010186677; 0.000138897];
    % Accelerometer noise standard deviation [m/s^2]. Models random, fluctuating errors.
    parameters.sensors.acc.std  = [0.011817097; 0.012306208; 0.015949601];
    
    %%% Magnetometer sensor
    % Earth's magnetic field vector in the inertial frame (e.g., pointing towards magnetic north).
    % This is the 'true' vector that the magnetometer measures to get a heading reference.
    parameters.sensors.mag.m_I  = [0.7071 ,0 ,0.7071]';  
    % Magnetometer constant bias [uT or normalized units].
    parameters.sensors.mag.bias = [0; 0; 0];
    % Magnetometer noise standard deviation [uT or normalized units].
    parameters.sensors.mag.std  = [0.011817097; 0.012306208; 0.015949601];
    
    %%% Gyroscope sensor
    % Gyroscope constant bias [rad/s]. Represents a slow drift in the angular rate measurement.
    parameters.sensors.gyro.bias = [-0.0000541; -0.0001553; -0.00006458];
    % Gyroscope noise standard deviation (Angular Random Walk) [rad/s]. Models high-frequency noise.
    parameters.sensors.gyro.std  = [0.009; 0.008; 0.007];
    % Gyroscope low-pass filter time constant [s]. Models the sensor's response time.
    parameters.sensors.gyro.tau = 0.1;
    
    %%% Star Tracker sensor (high-precision attitude sensor)
    % Star Tracker constant bias [rad].
    parameters.sensors.star.bias = 1E-6 * [1; 1; 1];
    % Star Tracker noise standard deviation [rad].
    parameters.sensors.star.std  = 1E-3 * [1; 1; 1];
    % Flag to enable (1) or disable (0) the star tracker in the simulation.
    parameters.sensors.star.enable = 1;
    % Number of stars (or direction vectors) to be simulated by the star tracker.
    parameters.sensors.star.numberOfStars = 10;
    
    % --- AHRS (Attitude and Heading Reference System) Selection ---
    % The AHRS is the algorithm that fuses sensor data to estimate the attitude.
    % Selector for the AHRS algorithm: 'MADGWICK', 'EKF', or 'UKF'.
    parameters.ahrs.flag = 'EKF';
    % Flag: 1 enables feedback from the AHRS to the controller (realistic).
    %       0 uses the ideal, true state for feedback (for baseline comparison).
    parameters.ahrs.enable = 1;
    
    % MADGWICK filter parameters
    parameters.ahrs.madwick.beta = 1; % Filter gain.
    parameters.ahrs.madwick.bias_gain = 0.01; % Gain for gyro bias estimation.
    
    % EKF (Extended Kalman Filter) parameters
    % The EKF needs a model of the sensor noise (R matrix).
    % Flag: 1 sets the EKF's noise parameters equal to the sensor model's noise (above).
    %       0 allows you to specify different noise values, e.g., if you want to
    %       intentionally mismatch the model to test robustness.
    parameters.ahrs.ekf.equalModel = 0;
    % EKF accelerometer noise standard deviation (used in R matrix).
    parameters.ahrs.ekf.acc_std    = [0.011817097; 0.012306208; 0.015949601];
    % EKF magnetometer noise standard deviation (used in R matrix).
    parameters.ahrs.ekf.mag_std    = [0.011817097; 0.012306208; 0.015949601];
    % EKF gyroscope noise standard deviation (used in Q matrix - process noise).
    parameters.ahrs.ekf.gyro_std   = [0.009; 0.008; 0.007];
    % EKF star tracker noise standard deviation (used in R matrix).
    parameters.ahrs.ekf.star_std   = 1E-3 * [1; 1; 1];
    
    % UKF (Unscented Kalman Filter) parameters
    % These parameters control the generation and weighting of sigma points.
    parameters.ahrs.ukf.alpha = 1e-3; % Scaling parameter, spread of sigma points.
    parameters.ahrs.ukf.beta  = 2.0;  % Used to incorporate prior knowledge of the state distribution.
    parameters.ahrs.ukf.kappa = 0.0;  % Secondary scaling parameter.
    parameters.ahrs.ukf.gyro_std = [0.009; 0.008; 0.007];
    parameters.ahrs.ukf.acc_std  = [0.011817097; 0.012306208; 0.015949601];
    parameters.ahrs.ukf.mag_std  = [0.011817097; 0.012306208; 0.015949601];
    parameters.ahrs.ukf.star_std  = 1E-3 * [1; 1; 1];
    
    % --- Controller Setpoint Parameters ---
    % The 'setpoint' is the target state that the controller tries to achieve.
    % Desired target attitude quaternion [w, x, y, z]'.
    parameters.setPoint.qd = [0.9808; 0.1691; 0.0933; 0.0277];
    % Desired target angular velocity [rad/s].
    parameters.setPoint.Wd = [0; 0; 0]; % e.g., hold a steady attitude.
    
    % --- Controller Selection ---
    % Selector: 1 for a standard PD Feedback Controller, 2 for Boskovic's Adaptive Controller.
    parameters.controller.selector = 1;
    
    % --- Feedback Controller (PD) Parameters ---
    % Derivative gain matrix (acts on angular velocity error).
    parameters.controller.feedback.Keye = 0.2 * eye(3);
    % Proportional gain matrix (acts on quaternion error).
    parameters.controller.feedback.Peye = 0.3 * eye(3);
    
    % --- Boskovic Adaptive Controller Parameters ---
    parameters.controller.boskController.delta = 0.5;
    parameters.controller.boskController.gamma = 0.001; % Adaptation rate.
    parameters.controller.boskController.Umax  = 0.10;  % Maximum control torque magnitude [Nm].
    parameters.controller.boskController.k0    = 1;     % Initial gain value.
    
    % --- Disturbance Torque Parameters ---
    % These model external torques acting on the satellite (e.g., solar pressure, drag).
    % Constant disturbance torque components [Nm].
    parameters.disturbances.constant.Tdx = 0.001;
    parameters.disturbances.constant.Tdy = 0.001;
    parameters.disturbances.constant.Tdz = 0.001;
    % Flag to enable (2) or disable (1) constant disturbances.
    parameters.disturbances.constant.enable = 1;
    
    % Parameters for variable (sinusoidal) disturbances [Nm].
    parameters.disturbances.variable.A   = 0.01;   % Amplitude
    parameters.disturbances.variable.w1  = 1.00;   % Frequency for x-axis
    parameters.disturbances.variable.w2  = 2.00;   % Frequency for y-axis
    parameters.disturbances.variable.w3  = 3.00;   % Frequency for z-axis
    % Flag to enable (2) or disable (1) variable disturbances.
    parameters.disturbances.variable.enable = 1;
    
    % --- Reaction Wheels (Actuators) Configuration ---
    % Number of reaction wheels in the assembly.
    parameters.rw.number = 3;
    % Maximum momentum storage per wheel [Nms]. A physical limitation.
    parameters.rw.maxH   = 1;
    % Mounting angles for wheels (beta1 and beta2 are typically used for cant angles).
    parameters.rw.beta1  = 0;
    parameters.rw.beta2  = 0;
    % Distribution matrix 'W'. Maps the 3-axis desired body torque to the
    % individual torque commands for each wheel based on their mounting orientation.
    % For a standard 3-orthogonal-wheel setup, this is the identity matrix.
    parameters.rw.W      = eye(3);
    
    % --- Reaction Wheel DC Motor Model ---
    % These parameters define the electrical and mechanical properties of the motors.
    parameters.rw.motor.kt   = 0.0254;   % Torque constant [Nm/A].
    parameters.rw.motor.Jrw  = 3.475E-5; % Rotor inertia [kg*m^2].
    parameters.rw.motor.b    = 5.58E-4;  % Viscous friction coefficient [Nms].
    parameters.rw.motor.c    = 0.0000;   % Coulomb friction (static friction) [Nm].
    parameters.rw.motor.L    = 0.000403; % Armature inductance [H].
    parameters.rw.motor.R    = 1.16;     % Armature resistance [Ohms].
    parameters.rw.motor.ke   = 0.0254;   % Back-EMF (voltage) constant [Vs/rad].
    
    % --- Reaction Wheel DC motor initial conditions ---
    parameters.rw.motor.init.w_rw = 0;      % Initial wheel speed [rad/s].
    parameters.rw.motor.init.current = 0; % Initial motor current [A].
    
    % --- Reaction Wheel PID Speed Controller Gains ---
    % This internal PID controller's job is to ensure the wheel's actual speed
    % matches the commanded speed from the ADCS control law.
    parameters.rw.motor.pid.kp = 0.504519552012394;  % Proportional gain
    parameters.rw.motor.pid.ki = 28.7345189644697;   % Integral gain
    parameters.rw.motor.pid.kd = -0.000554999924740984; % Derivative gain
    parameters.rw.motor.pid.Nu = 289.873462888377;  % Derivative filter coefficient.
end