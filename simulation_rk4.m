function [x, u, T_winf_nosat, x_est, indicators, sensors, error_flag] = simulation_rk4(app,disturbances,simParameters,time)
% SIMULATION_RK4_WITH_SAMPLING Simulates CubeSat attitude using a 4th-order Runge-Kutta integrator and EKF, including sensor sampling times.
%
% This modified version of the program simulates the attitude dynamics of a
% CubeSat, allowing sensors to operate at different sampling rates. The
% gyroscope, used for EKF prediction, can operate at a high frequency,
% while vector sensors (accelerometer, magnetometer, star tracker), used
% for correction, are sampled at a lower frequency.
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
%                   - .sensors.attitude.Ts: Sampling time for attitude sensors
%                     (accelerometer, magnetometer, etc.) (s).
%   time          - Structure defining simulation time parameters.
%
% Outputs:
%   x             - Actual CubeSat state over time.
%   u             - Control torque applied by the actuator.
%   x_est         - EKF estimated CubeSat state over time.
%   indicators    - Structure containing performance metrics.
%   sensors       - Structure containing sensor data. Measurements not
%                   available at a time step are marked as NaN.
%   error_flag    - Error flag: 1 if NaN was detected.
%
% Detailed Description of Changes:
%   - **Sampling Times**: Two key parameters are introduced in `simParameters.sensors`:
%     `Ts_gyro` for the gyroscope and `Ts_attitude` for the vector sensors.
%   - **Sampling Logic**: Inside the main loop, the `mod` operator is used to
%     determine if it's time to take a new sensor sample.
%   - **Data Hold**: If a new sample is not taken, the value of the
%     previous measurement is held and used for calculations (Zero-Order Hold).
%     This is especially important for the gyroscope, which feeds the
%     EKF prediction at every step.
%   - **EKF Predict-Update**: The EKF **prediction** step is executed at every
%     simulation step (`dt`) using the most recent gyroscope reading.
%     The **correction/update** step is only executed when a new
%     batch of attitude sensor measurements is available. In the intermediate
%     steps, the state estimate is simply the result of the prediction.
% ----------------------------------------------------------------------------------
%%
simParameters.sensors.gyro.Ts = 5/1000;      % Gyro sampling time [s] -> 100 Hz
simParameters.sensors.attitude.Ts = 10/1000; % Attitude sensors sampling time [s] -> 1 Hz
simParameters.control.Ts = 1/10;

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
if app.controller_popupmenu.ValueIndex == 1
    %%% Feedback Controller gains
    P = simParameters.feedback.Peye; 
    K = simParameters.feedback.Peye;
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
m_I      = [0,1,0]';  % Magnetometer reference direction (in inertial frame)
std_mag  = simParameters.sensors.mag.std;  % Standard deviation
bias_mag = simParameters.sensors.mag.bias; % Bias of magnetometer
%%% Gyroscope sensor
bias_gyro = simParameters.sensors.gyro.bias;  % Gyroscope bias
std_gyro  = simParameters.sensors.gyro.std;   % Standard deviation

% Retrieve sampling times from the parameters structure.
Ts_gyro     = simParameters.sensors.gyro.Ts;     % Sampling time for gyro
Ts_attitude = simParameters.sensors.attitude.Ts; % Sampling time for correction
Ts_control  = simParameters.control.Ts;          % Sampling time for control 

% Calculate in how many simulation steps (dt) each sampling occurs.
% round() is used to find the nearest integer number of steps.
steps_gyro     = round(Ts_gyro / dt);
steps_attitude = round(Ts_attitude / dt);
steps_control  = round(Ts_control / dt);
% Ensure at least one step
if steps_gyro < 1, steps_gyro = 1; end
if steps_attitude < 1, steps_attitude = 1; end
if steps_control < 1, steps_control = 1; end

%%% Star Sensor
if simParameters.sensors.star.enable == 1  % Check if star sensor is enabled
    number_of_stars = simParameters.sensors.star.numberOfStars;  % Number of stars to be used by the sensor
    % Generate random azimuth and elevation angles for stars
    theta = 2*pi*rand(1, number_of_stars);       % Azimuth angles
    phi = acos(2*rand(1, number_of_stars) - 1);  % Elevation angles     
    % Convert from spherical to Cartesian coordinates
    x_I = sin(phi) .* cos(theta);  % X component of star vectors
    y_I = sin(phi) .* sin(theta);  % Y component of star vectors
    z_I = cos(phi);                % Z component of star vectors
    % Create a 3xN matrix of unit vectors representing the stars' positions in space
    stars_I = [x_I; y_I; z_I];  % Matrix of unit vectors for stars
    bias_star = simParameters.sensors.star.bias;  % Bias
    std_star  = simParameters.sensors.star.std;   % Standard deviation
else
    number_of_stars = 0;  % If star sensor is not enabled, set the number of stars to 0
    stars_I = [];
end
%%% EKF gains
%%% Check if the sensor model std is equal to the EKF model
if simParameters.ekf.equalModel
    if simParameters.sensors.star.enable == 1
        std_star_ekf = std_star;
    end
    std_acc_ekf  = std_acc;
    std_mag_ekf  = std_mag;
    std_gyro_ekf = std_gyro;
else
    if simParameters.sensors.star.enable == 1
        std_star_ekf = simParameters.ekf.star.std;
    end
    std_acc_ekf  = simParameters.ekf.acc.std;
    std_mag_ekf  = simParameters.ekf.mag.std;
    std_gyro_ekf = simParameters.ekf.gyro.std;
end
%%% Previous calculus to define Q
Q_gyro = diag(std_gyro_ekf.^2);
%%% Initial state covariance matrix P (7x7 identity matrix)
P_cov_ant = eye(7);  
%%% Define measurements covariance matrix R
%%% If star sensor is enabled, include the star sensor noise in the covariance matrix
if simParameters.sensors.star.enable == 1
    star = diag(std_star_ekf.^2);
    star_large = kron(eye(number_of_stars), star);
    R_k = blkdiag(diag(std_acc_ekf.^2), diag(std_mag_ekf.^2), star_large);
else
    R_k = blkdiag(diag(std_acc_ekf.^2), diag(std_mag_ekf.^2));
end
%%% Actuator parameters
number_of_rw = simParameters.rw.number;
W  = simParameters.rw.W;
%% 2. Simulation containers
%%% System states containers
x = NaN(7,n); u = NaN(3,n); dq = NaN(4,n-1); T_winf_nosat = NaN(number_of_rw,n);
%%% Sensor containers
g_B = NaN(3,n-1); m_B = NaN(3,n-1); stars_B = NaN(3*number_of_stars,n-1);
omega_meas = NaN(3,n-1);
%%% Kalman filter containers
x_est = NaN(7,n); y_est = NaN(6+3*number_of_stars,n-3); 
%%% Performance index containers
o = NaN(1,n-1); 
%% 3. Initial conditions
x(:,1)=[simParameters.initialValues.q0; simParameters.initialValues.Wo];
u(:,1) = zeros(3,1);  x_est(:,1) = [1, zeros(1,6)]'; 
T_winf_nosat(:,1) = zeros(number_of_rw, 1);
%Initialize last known attitude sensor measurements
R_init = my_quat2rot(x(1:4, 1));
last_gB_meas = accelerometer_readings(R_init, g_I, std_acc, bias_acc); % Initial noise-free reading
last_mB_meas = magnetometer_measurements(R_init, m_I, std_mag, bias_mag); % Initial noise-free reading
if simParameters.sensors.star.enable == 1 
    last_starSens_meas = star_sensors_measurement(R_init,stars_I,std_star,bias_star,number_of_stars);
end
% Initialize containers for the last known measurement (Zero-Order Hold).
last_omega_meas = simParameters.initialValues.Wo; % Initialize with the initial real angular velocity.
last_u = zeros(3,1); last_T_winf_nosat = zeros(number_of_rw, 1);

%% 4. Previous calculations to optimize the simulation
%%%Calculate Wd_dot
Wd_dot = diff(wd')'./diff(t);
%%%Create wait bar
hWaitbar = waitbar(0, 'Progress: 0%','Name', 'Attitude simulation progress...');
%% 5. Initialize simulation
for i = 1:n-1
    %% 5.1. Sensor models (add noise, bias and SAMPLING)
    
    % Gyroscope model with sampling.
    % Check if the current step 'i' is a sampling instant for the gyroscope.
    if mod(i-1, steps_gyro) == 0
        W_gyro = std_gyro .* randn(3, 1);                  % Random noise
        current_omega_meas = x(5:7, i) + W_gyro + bias_gyro; % Gyroscope measurement
        last_omega_meas = current_omega_meas; % Update the last known measurement.
    else
        current_omega_meas = last_omega_meas; % Hold the previous measurement.
    end
    omega_meas(:, i) = current_omega_meas; % Save to history.

    %% 5.2. EKF - Extended Kalman Filter
    
    % --- EKF Initial Guess (TRIAD at first step) ---
    if i == 1
        % Force an initial measurement for the TRIAD algorithm.
        current_gB_meas = last_gB_meas;
        current_mB_meas = last_mB_meas;
        if simParameters.sensors.star.enable == 1 
            current_starSens_meas = last_starSens_meas;
        end
        q_est_init = triad_algorithm(g_I, m_I, current_gB_meas, current_mB_meas);
        x_est(1:4, i) = q_est_init;
    end
    
    % --- EKF Prediction Step (always runs) ---
    % The prediction step runs at each time step 'dt' using the
    % most recent gyroscope measurement (which may be a held measurement).
    A = [eye(4), -0.5 *dt* xi_matrix(x_est(1:4, i)); zeros(3, 4), eye(3)];
    B = 1/2 * dt * [xi_matrix(x_est(1:4, i)); zeros(3, 3)];
    x_pred = A * x_est(:, i) + B * omega_meas(:, i);
    
    Q = dt*[0.25*xi_matrix(x_est(1:4,i))*Q_gyro*xi_matrix(x_est(1:4,i))', zeros(4,3);
              zeros(3,4), Q_gyro];
    P_pred = A * P_cov_ant * A' + Q; 

    % The EKF correction step only runs when new attitude measurements are available.
    if mod(i-1, steps_attitude) == 0
        % Generate new attitude sensor measurements
        R = my_quat2rot(x(1:4, i));
        % Accelerometer readings
        current_gB_meas = accelerometer_readings(R, g_I, std_acc, bias_acc);
        last_gB_meas = current_gB_meas; % Update the last known measurement.

        % Magnetometer readings
        current_mB_meas = magnetometer_measurements(R, m_I, std_mag, bias_mag);
        last_mB_meas = current_mB_meas; % Update the last known measurement.
        
        % Star sensor model
        if simParameters.sensors.star.enable == 1
            current_starSens_meas = star_sensors_measurement(R,stars_I,std_star,bias_star,number_of_stars);
            last_starSens_meas = current_starSens_meas;
            y = [last_gB_meas; last_mB_meas; last_starSens_meas];
        else
            y = [last_gB_meas; last_mB_meas];
        end
        
        % --- EKF Correction Step ---
        C = measurement_jacobian(x_pred(1:4), g_I, m_I, stars_I); % Jacobian evaluated at the predicted state
        y_est_pred = C * [x_pred(1:4); zeros(3,1)]; % Predicted measurement (only depends on the quaternion)
        y_est(:,i) = y_est_pred;

        k_K = P_pred * C' / (C * P_pred * C' + R_k);  % Kalman gain
        x_est(:, i + 1) = x_pred + k_K * (y - y_est(:,i)); % Correct state
        P_cov = (eye(7) - k_K * C) * P_pred; % Correct covariance
    else
        % --- No new measurement: propagate prediction ---
        % If there are no new attitude measurements, the estimated state for the
        % next step is simply the predicted state.
        x_est(:, i + 1) = x_pred;
        P_cov = P_pred;
    end
    
    % --- Normalize quaternion and update covariance ---
    x_est(1:4, i + 1) = x_est(1:4, i + 1) / norm(x_est(1:4, i + 1));
    if x_est(1, i + 1) < 0
         x_est(1:4, i + 1) = -x_est(1:4, i + 1);
    end
    P_cov_ant = P_cov;
    
    %%% Store sensors measurements
    g_B(:,i) = current_gB_meas;
    m_B(:,i) = current_mB_meas;
    if simParameters.sensors.star.enable == 1
        stars_B(:,i) = current_starSens_meas;
    end

    %% 5.3. Runge-Kutta 4th order integration for state update
    T_u = W*T_winf_nosat(:, i);
    g1 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i));
    g2 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g1);
    g3 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g2);
    g4 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g3);
    x(:, i + 1) = x(:, i) + (1 / 6) * (g1 + 2 * g2 + 2 * g3 + g4);
    
    % Calculate quaternion error
    feed_est = [x_est(1:4,i); omega_meas(:,i)]; % Use the most recent omega measurement for control
    dq(:, i) = Error_quaternio(qd(:, i), feed_est(1:4));
    
    % Control law and computational cost estimation
    tic
    if app.controller_popupmenu.ValueIndex == 1
        u(:, i + 1) = ControlFeedback_rw(I, feed_est, dq(:, i), wd(:, i), Wd_dot(:, i), P, K); 
    else
        k_dot = Gain_estimator_bosk(feed_est(5:7), wd(:, i), dq(:, i), delta, gamma, k, Umax);
        k = k_ant + dt / 6 * (k_dot_ant + 2 * (k_dot_ant + k_dot) + k_dot);
        u(:, i + 1) = Boskovic_control(feed_est(5:7), wd(:, i), dq(:, i), delta, k, Umax);
    end
    o(i) = toc;
    
    % Allocator
    if number_of_rw == 3
        T_winf_nosat(:, i + 1) = W\u(:, i + 1);
    else
        T_w_L2norm = allocator_L2norm(W, u(:, i + 1));
        T_winf_nosat(:, i + 1) = allocator_LinfNorm(T_w_L2norm); 
    end
    % Check for NaN errors
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
end
sensors.meas = [omega_meas; g_B; m_B; stars_B];
sensors.est = y_est;
%%%Close wait bar
close(hWaitbar);
%% 6. Update performance indices
    if error_flag == 0
        %%%EULERINT calculation
        ang_and_axis = quat2axang(dq'); 
        eulerang = ang_and_axis(:,4);
        indicators.eulerInt = cumtrapz(t(1:end-1),eulerang); 
        %%%ASCCT Calculation
        ascct_dot=vecnorm(u,2,1).^2;
        indicators.ascct = cumtrapz(t,ascct_dot); 
        %%%Settlement time calculation
        tol = 5/100; % 5% of the set point
        indicators.ts = calculateSettlementTime(180/pi*quat2eul(dq'), t, tol);
        %%%Computation time
        indicators.o = o;
    end
end
%% 7. Program Functions
function Xi = xi_matrix(q)
% xi_matrix calculates a 4x3 matrix from a 4-element vector q.
% Developed by bespi123
    Xi = [ -q(2), -q(3), -q(4);
            q(1), -q(4),  q(3);
            q(4),  q(1), -q(2);
           -q(3),  q(2),  q(1) ];
end
function R = my_quat2rot(q)
% my_quat2rot Converts a quaternion to a 3x3 rotation matrix.
% Developed by bespi123
    q = q / norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end
function C = measurement_jacobian(q, g_I, m_I, star_I)
% This function calculates the measurement Jacobian matrix for an EKF.
% Developed by bespi123
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    gx = g_I(1); gy = g_I(2); gz = g_I(3);
    rx = m_I(1); ry = m_I(2); rz = m_I(3);
    C_g = 2*[gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, -gx*q2 + gy*q1 - gz*q0, -gx*q3 + gy*q0 + gz*q1, 0, 0, 0;
            -gx*q3 + gy*q0 + gz*q1,  gx*q2 - gy*q1 + gz*q0,  gx*q1 + gy*q2 + gz*q3, -gx*q0 - gy*q3 + gz*q2, 0, 0, 0;
             gx*q2 - gy*q1 + gz*q0,  gx*q3 - gy*q0 - gz*q1,  gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, 0, 0, 0];
    C_m = 2*[rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, -rx*q2 + ry*q1 - rz*q0, -rx*q3 + ry*q0 + rz*q1, 0, 0, 0;
            -rx*q3 + ry*q0 + rz*q1,  rx*q2 - ry*q1 + rz*q0,  rx*q1 + ry*q2 + rz*q3, -rx*q0 - ry*q3 + rz*q2, 0, 0, 0;
             rx*q2 - ry*q1 + rz*q0,  rx*q3 - ry*q0 - rz*q1,  rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, 0, 0, 0];
    stars_num = size(star_I,2);
    C_star = NaN(3*stars_num,7);
    for i=1:stars_num
        star_x = star_I(1,i); star_y = star_I(2,i); star_z = star_I(3,i);
        index1 = ((i-1)*3+1); index2 = (index1+3)-1;
        C_star(index1:index2,:) = 2*[star_x*q0 + star_y*q3 - star_z*q2,  star_x*q1 + star_y*q2 + star_z*q3, -star_x*q2 + star_y*q1 - star_z*q0, -star_x*q3 + star_y*q0 + star_z*q1, 0, 0, 0;
            -star_x*q3 + star_y*q0 + star_z*q1,  star_x*q2 - star_y*q1 + star_z*q0,  star_x*q1 + star_y*q2 + star_z*q3, -star_x*q0 - star_y*q3 + star_z*q2, 0, 0, 0;
             star_x*q2 - star_y*q1 + star_z*q0,  star_x*q3 - star_y*q0 - star_z*q1,  star_x*q0 + star_y*q3 - star_z*q2,  star_x*q1 + star_y*q2 + star_z*q3, 0, 0, 0];
    end
    C = [C_g; C_m; C_star];
end
function q_triad = triad_algorithm(r1,r2,b1,b2)
% triad_algorithm Calculates the attitude matrix using the TRIAD algorithm.
% Developed by bespi123
    r1 = reshape(r1, [], 1); r2 = reshape(r2, [], 1); 
    b1 = reshape(b1, [], 1); b2 = reshape(b2, [], 1);
    v2 = cross(r1,r2) / norm(cross(r1,r2));
    w2 = cross(b1,b2) / norm(cross(b1,b2));
    A_triad = b1*r1' + (cross(b1,w2))*(cross(r1,v2))' + w2*v2';
    q_triad = dcm2quat(A_triad)';
end
function g_B = accelerometer_readings(R,g_I,std_acc,bias_acc)
    % Accelerometer readings
    W_acc = std_acc .* randn(3, 1);
    g_B = R' * g_I + W_acc + bias_acc;
    g_B = g_B/norm(g_B);
end

function m_B = magnetometer_measurements(R,m_I,std_mag,bias_mag)
    W_mag = std_mag .* randn(3, 1);
    m_B = R' * m_I + W_mag + bias_mag;
    m_B = m_B/norm(m_B);
end

function stars_B = star_sensors_measurement(R,stars_I,std_star,bias_star,number_of_stars)
    W_star = std_star .* rand(3,number_of_stars);
    stars_B_temp = R' * stars_I + W_star + bias_star;
    stars_B_temp = stars_B_temp ./ vecnorm(stars_B_temp, 2);
    stars_B = stars_B_temp(:);
end