function [x, u, T_winf_nosat, x_est, indicators, sensors, actuators, error_flag] = simulation_rk4(~,disturbances,simParameters,time)
% SIMULATION_RK4 Simulates CubeSat attitude with sampling for sensors and control loop.
%
% This version simulates the attitude dynamics of a CubeSat, allowing the
% sensors and the control loop to operate at different sampling rates, which
% is a more realistic model of a digital control system.
%
% Author: bespi123
% Creation Date: [Creation Date, 2023-10-27]
% Last Modified: [Last Modification Date, 2024-07-08]
%
% Inputs:
%   app           - UI application object.
%   disturbances  - Vector of external disturbance torques.
%   simParameters - Structure containing simulation parameters, including:
%                   - .sensors.gyro.Ts: Gyroscope sampling time (s).
%                   - .sensors.attitude.Ts: Sampling time for attitude sensors (s).
%                   - .control.Ts: Sampling time for the control loop (s).
%   time          - Structure defining simulation time parameters.
%
% Outputs:
%   x             - Actual CubeSat state over time.
%   u             - Control torque applied (updated according to Ts_control).
%   x_est         - EKF estimated CubeSat state over time.
%   indicators    - Structure containing performance metrics.
%   sensors       - Structure containing sensor data.
%   error_flag    - Error flag: 1 if NaN was detected.
%
% Detailed Description of Changes:
%   - **Control Sampling**: A `simParameters.control.Ts` is introduced to define
%     how often the control law is computed.
%   - **Zero-Order Hold (ZOH)**: The calculated control torque is held constant
%     between control loop updates, simulating a real digital system's behavior.
%   - **EKF Predict-Update**: The EKF logic remains: the **prediction** step
%     runs at the high frequency of the gyroscope, while the **correction** step
%     runs at the lower frequency of the attitude sensors. This is the standard
%     and most effective implementation.
% ----------------------------------------------------------------------------------
%% Imports
%%% Import adsim utilities
import adcsim.utils.*
%%% Import environment to use python libraries
%pyenv("Version", "D:\git\Satellite_ADCS_GUIDE\.venv\Scripts\python.exe"); 

%% 1. Parameter Recovery and Initialization
%%% Time parameters assuming uniform time steps
t = time.t; n = time.n; dt = t(2) - t(1);
%%% Desired states
qd = repmat(simParameters.setPoint.qd',n, 1)';    
wd = repmat(simParameters.setPoint.Wd',n, 1)';
%%% Disturbances and Inertia
Td = disturbances; I = simParameters.initialValues.I;
%%% Error flag
error_flag = 0;
% Controller gains and initial values
if simParameters.controller.selector == 1
    %%% Feedback Controller gains
    P = simParameters.feedback.Peye; 
    K = simParameters.feedback.Keye;
else
    %%% Boskovic Controller gains
    delta = simParameters.boskController.delta;
    gamma = simParameters.boskController.gamma;
    k     = simParameters.boskController.k0;
    Umax  = simParameters.boskController.Umax;
    %%% Previous gains
    k_ant = k; k_dot_ant = 0;
end
% Sensor model parameters (inertial frame references)
%%% Accelerometer sensor
g_I      = [0,0,1]';  % Gravitational acceleration vector in inertial frame 
std_acc  = simParameters.sensors.acc.std;  % Standard deviation
bias_acc = simParameters.sensors.acc.bias; % Bias of accelerometer
%%% Magnetometer sensor
m_I      = [0.7071 ,0 ,0.7071]';  % Magnetometer reference direction (in inertial frame)
std_mag  = simParameters.sensors.mag.std;  % Standard deviation
bias_mag = simParameters.sensors.mag.bias; % Bias of magnetometer
%%% Gyroscope sensor
bias_gyro = simParameters.sensors.gyro.bias;  % Gyroscope bias
std_gyro  = simParameters.sensors.gyro.std;   % Standard deviation

%%% Gyro filter time constant
simParameters.sensors.gyro.tau = 0.1;

% Retrieve sampling times from the parameters structure.
Ts_gyro = simParameters.sensors.gyro.Ts;
Ts_attitude = simParameters.sensors.attitude.Ts;
Ts_control = simParameters.control.Ts; % Sampling time for the control loop

% Calculate in how many simulation steps (dt) each sampling occurs.
steps_gyro = round(Ts_gyro / dt);
steps_attitude = round(Ts_attitude / dt);
steps_control = round(Ts_control / dt); % Steps for the control loop
if steps_gyro < 1, steps_gyro = 1; end % Ensure at least one step
if steps_attitude < 1, steps_attitude = 1; end
if steps_control < 1, steps_control = 1; end

%%% Star Sensor
if simParameters.sensors.star.enable == 1
    number_of_stars = simParameters.sensors.star.numberOfStars;
%    % Generate random star vectors in the inertial frame
%    theta = 2*pi*rand(1, number_of_stars);
%    phi = acos(2*rand(1, number_of_stars) - 1);
%    x_I = sin(phi) .* cos(theta); y_I = sin(phi) .* sin(theta); z_I = cos(phi);
%    stars_I = [x_I; y_I; z_I];
    bias_star = simParameters.sensors.star.bias;
    std_star  = simParameters.sensors.star.std;
%else
%    number_of_stars = 0;
%    stars_I = [];
end

%%% EKF gains
if simParameters.ekf.equalModel
    if simParameters.sensors.star.enable == 1, std_star_ekf = std_star; end
    std_acc_ekf  = std_acc; std_mag_ekf  = std_mag; std_gyro_ekf = std_gyro;
else
    if simParameters.sensors.star.enable == 1, std_star_ekf = simParameters.ekf.star.std; end
    std_acc_ekf  = simParameters.ekf.acc.std; std_mag_ekf  = simParameters.ekf.mag.std; std_gyro_ekf = simParameters.ekf.gyro.std;
end
%%% Process noise covariance for gyro
Q_gyro = diag(std_gyro_ekf.^2);
%%% Initial state covariance matrix P
P_cov_ant = eye(7);  
%%% Define measurements covariance matrix R
if simParameters.sensors.star.enable == 1
    star = diag(std_star_ekf.^2);
    star_large = kron(eye(number_of_stars), star);
    R_k = blkdiag(diag(std_acc_ekf.^2), diag(std_mag_ekf.^2), star_large);
else
    number_of_stars = 0;
    R_k = blkdiag(diag(std_acc_ekf.^2), diag(std_mag_ekf.^2));
end

%%% Actuator parameters
number_of_rw = simParameters.rw.number;
W  = simParameters.rw.W;

%%% Brushless models Motor Parameters
motor.kt   = simParameters.rw.motor.kt;  %kt [N*m/A] Torque constant
motor.Jrw  = simParameters.rw.motor.Jrw; %J  [kg*m^2] Rotor inertia
motor.b    = simParameters.rw.motor.b;   %B  [N*m*s] Viscous friction
motor.c    = simParameters.rw.motor.c;   %kc [N*m] Coulomb friction
motor.L    = simParameters.rw.motor.L;   %H  [H] Inductance
motor.R    = simParameters.rw.motor.R;   %R  [Ω] Resistance
motor.ke   = simParameters.rw.motor.ke;  %ke [V*s/rad] Back-EMF constant

%%Initial conditions
init.w_rw = 0;
init.current = 0; 

%% 2. Simulation containers
x = NaN(7,n); u = NaN(3,n); dq = NaN(4,n-1); T_winf_nosat = NaN(number_of_rw,n);
g_B = NaN(3,n-1); m_B = NaN(3,n-1); stars_B = NaN(3*number_of_stars,n-1);
omega_meas = NaN(3,n-1); omega_meas_filtered = NaN(3,n-1);
x_est = NaN(7,n-1); y_est = NaN(6+3*number_of_stars,n-3); 
o = NaN(1,n-1); 

x_rw = NaN(2*number_of_rw,n);
u_rw = NaN(number_of_rw,n);
w_cmd = NaN(number_of_rw,n);
torque_real = NaN(number_of_rw,n);

%% 3. Initial conditions
% Set initial true state, control, and estimated state
x(:,1)=[simParameters.initialValues.q0; simParameters.initialValues.Wo];
u(:,1) = zeros(3,1);  x_est(:,1) = [1, zeros(1,6)]'; 
T_winf_nosat(:,1) = zeros(number_of_rw, 1);

x_rw(:,1) = repmat([init.w_rw, init.current],1,number_of_rw)';
u_rw(:,1) = zeros(number_of_rw,1);
w_cmd_ant = zeros(number_of_rw,1);

% Initialize containers for the last known measurement/calculation (Zero-Order Hold).
last_u = zeros(number_of_rw,1);
last_T_winf_nosat = zeros(number_of_rw, 1);

% Initialize simulation objects
mySatellite = adcsim.satellite.Satellite(simParameters.initialValues, I);
reactionWheels = adcsim.actuators.ReactionWheelAssembly(simParameters.rw);
myIMU =  adcsim.sensors.IMU(simParameters.sensors, g_I, m_I);
myStarSensor = adcsim.sensors.StarTracker(simParameters.sensors.star);
myEKF = adcsim.AHRS_algorithms.AttitudeEKF(simParameters.ekf, dt, myIMU, myStarSensor);
madwick_params.beta = 1;
madwick_params.bias_gain = 0.01;
myMadwick = adcsim.AHRS_algorithms.MadgwickAHRS(madwick_params, dt, myIMU, myStarSensor);

ukf_params.alpha = 1e-3;
ukf_params.beta  = 2.0;
ukf_params.kappa = 0.0;
ukf_params.gyro.std = std_gyro;
ukf_params.acc.std  = std_acc;
ukf_params.mag.std  = std_mag;
ukf_params.star.std  = std_star;

myUKF = adcsim.AHRS_algorithms.AttitudeUKF(ukf_params, dt, myIMU, myStarSensor);

% Get the initial rotation matrix from the true state.
R_init = my_quat2rot(x(1:4, 1));

% Simulate an initial, sensor reading. 
g_B(:,1) = myIMU.getAccelerometerReading(R_init); 
m_B(:,1) = myIMU.getMagnetometerReading(R_init); 
myEKF.initializeWithTriad(g_B(:,1), m_B(:,1));

if simParameters.sensors.star.enable == 1
    stars_B(:,1) = myStarSensor.getReading(R_init);
end

%simParameters.ahrs.flag = 'MADGWICK';
%simParameters.ahrs.flag = 'EKF';
simParameters.ahrs.flag = 'UKF';

%% 4. Previous calculations to optimize the simulation
Wd_dot = diff(wd')'./diff(t);
hWaitbar = waitbar(0, 'Progress: 0%','Name', 'Attitude simulation progress..');
%% 5. Initialize simulation
for i = 1:n-1
    %% 5.1. Sensor models (add noise, bias and SAMPLING)
    % Gyroscope model with sampling.
    if mod(i-1, steps_gyro) == 0 || i==1
        current_omega_meas  = myIMU.getGyroscopeReading(x(5:7, i));
        filtered_omega_meas = myIMU.filterGyroscopeReading(current_omega_meas, Ts_gyro);
    end
    omega_meas(:, i) = current_omega_meas;          % Save to history.
    omega_meas_filtered(:,i) = filtered_omega_meas; % Save filtered data to history
    
    %% 5.2. EKF - Extended Kalman Filter
    if strcmp(simParameters.ahrs.flag, 'EKF')
        % --- EKF Prediction Step (always runs at high frequency) ---
        % The prediction step runs at each time step 'dt' using the
        % most recent gyroscope measurement (which may be a held measurement).
        myEKF.predict(filtered_omega_meas);

        %%% --- EKF Correction Step (runs only when attitude measurements are available) ---
        if mod(i-1, steps_attitude) == 0
            % Generate new attitude sensor measurements using helper functions
            R = my_quat2rot(x(1:4, i));
            g_B(:,i) = myIMU.getAccelerometerReading(R);
            m_B(:,i) = myIMU.getMagnetometerReading(R);

            if simParameters.sensors.star.enable == 1 
                stars_B(:,i)  = myStarSensor.getReading(R);
                y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
            else
                y = [g_B(:,i); m_B(:,i)];
            end

            % Perform the correction
            myEKF.correct(y); 
        end
        x_est(:, i + 1) = myEKF.x_est';
    elseif strcmp(simParameters.ahrs.flag, 'MADGWICK')
        %% Madwick Algorithm
        % Generate new attitude sensor measurements using helper functions
        myMadwick = myMadwick.Predict(filtered_omega_meas);
        if mod(i-1, steps_attitude) == 0
            R = my_quat2rot(x(1:4, i));
            g_B(:,i) = myIMU.getAccelerometerReading(R);
            m_B(:,i) = myIMU.getMagnetometerReading(R);
            if simParameters.sensors.star.enable == 1 
                stars_B(:,i)  = myStarSensor.getReading(R);
                y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
            else
                y = [g_B(:,i); m_B(:,i)];
            end
            myMadwick = myMadwick.Update(filtered_omega_meas, y);
        end
       
        x_est(1:4, i + 1) = myMadwick.x_est;
        x_est(5:7, i + 1) = myMadwick.gyro_bias;
    else
       % Generate new attitude sensor measurements using helper functions
        myUKF.Predict(filtered_omega_meas);
        if mod(i-1, steps_attitude) == 0
            R = my_quat2rot(x(1:4, i));
            g_B(:,i) = myIMU.getAccelerometerReading(R);
            m_B(:,i) = myIMU.getMagnetometerReading(R);
            if simParameters.sensors.star.enable == 1 
                stars_B(:,i)  = myStarSensor.getReading(R);
                y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
            else
                y = [g_B(:,i); m_B(:,i)];
            end
            myUKF.Correct(y);
        end


       x_est(1:4, i + 1) = myUKF.x_est;
       x_est(5:7, i + 1) = zeros(3,1);
    end
    
    %% 5.3. Control Law (runs at its own sampling rate Ts_control)
    if mod(i-1, steps_control) == 0
        % --- Calculate new control command ---
        if simParameters.ekf.enable == 1
            %%% Use sensors model as feedback signal
            feed_est = [x_est(1:4,i); omega_meas_filtered(:,i)-x_est(5:7,i)];
        else
            %%% Use ideal signals as feedback
            feed_est = x(:, i);
        end
        
        dq(:, i) = Error_quaternio(qd(:, i), feed_est(1:4));
        
        tic
        if simParameters.controller.selector == 1
            u_new = adcsim.controllers.ControlFeedback_rw(I, feed_est, dq(:, i), wd(:, i), Wd_dot(:, i), P, K); 
        else
            k_dot = adcsim.controllers.Gain_estimator_bosk(feed_est(5:7), wd(:, i), dq(:, i), delta, gamma, k, Umax);
            %k = k_ant + dt / 6 * (k_dot_ant + 2 * (k_dot_ant + k_dot) + k_dot);
            k = k_ant + Ts_control / 6 * (k_dot_ant + 2 * (k_dot_ant + k_dot) + k_dot);
            u_new = adcsim.controllers.Boskovic_control(feed_est(5:7), wd(:, i), dq(:, i), delta, k, Umax);
        end
        o(i) = toc;
        
        % --- Allocator ---
        if number_of_rw == 3
            T_winf_new = W\u_new;
        else
            T_w_L2norm = allocator_L2norm(W, u_new);
            T_winf_new = allocator_LinfNorm(T_w_L2norm); 
        end
        
        % --- Update last known values (for ZOH) ---
        last_u = u_new;
        last_T_winf_nosat = T_winf_new;

        % --- Get rw angular rate command
        last_w_rw_cmd = w_cmd_ant + (last_T_winf_nosat / motor.Jrw) * Ts_control;

    else
        % --- Hold previous control command (ZOH) ---
        o(i) = NaN; % No control computation cost
    end
    
    % Log the control signal for every step (it will be piecewise constant)
    u(:, i + 1) = last_u;
    T_winf_nosat(:, i + 1) = last_T_winf_nosat;
    w_cmd(:, i) = last_w_rw_cmd;


    %% Actuator model
    [reactionWheels, torque_real_vector] = reactionWheels.update(w_cmd(:, i), dt);
    x_rw(:, i + 1) = reactionWheels.States(:);
    torque_real(:,i) = torque_real_vector';
    T_u  = W*torque_real(:,i);
    mySatellite = mySatellite.updateState(T_u, Td(:, i), dt);
    x(:, i + 1) = mySatellite.State;

    if isnan(x(:, i + 1))
        error_flag = 1;
        break;
    end
    %% Progress bar update
    if mod(i, round(n / 10)) == 0
        progress = i / (n - 1);
        if isa(hWaitbar, 'handle') && isvalid(hWaitbar)
            waitbar(progress, hWaitbar, sprintf('Progress: %.1f%%', progress * 100));
        end
    end

    w_cmd_ant = w_cmd(:,i);
end
sensors.meas = [omega_meas; fillmissing([g_B; m_B; stars_B], 'previous',2)];
sensors.est = fillmissing(y_est, 'previous',2);
sensors.omega_filtered = omega_meas_filtered;

actuators.w_cmd = w_cmd;
actuators.x_rw  = x_rw;
actuators.torque_real = torque_real;
actuators.u_rw  = u_rw;

close(hWaitbar);
%% 6. Update performance indices
    if error_flag == 0
        ang_and_axis = quat2axang(fillmissing(dq, 'previous',2)'); 
        eulerang = ang_and_axis(:,4);
        indicators.eulerInt = cumtrapz(t(1:end-1),eulerang); 
        ascct_dot=vecnorm(u,2,1).^2;
        indicators.ascct = cumtrapz(t,ascct_dot); 
        tol = 5/100;
        indicators.ts = calculateSettlementTime(180/pi*quat2eul(dq'), t, tol);
        indicators.o = fillmissing(o, 'previous');
        
        ang_and_axis_real = quat2axang(x(1:4,:)');
        ang_and_axis_est = quat2axang(x_est(1:4,:)');
        error_vec = ang_and_axis_real(:,4) - ang_and_axis_est(:,4);
        indicators.RMSE = sqrt(mean(error_vec.^2));         
    else
        errordlg('Unstable System', 'ERROR');
        indicators = NaN;
    end
end
% % % %% 7. Program Functions
function R = my_quat2rot(q)
% my_quat2rot Converts a quaternion to a 3x3 rotation matrix.
% Developed by bespi123
    q = q / norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end