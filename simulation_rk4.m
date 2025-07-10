function [x, u, T_winf_nosat, x_est, indicators, sensors, error_flag] = simulation_rk4(~,disturbances,simParameters,time)
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
m_I      = [0,1,0]';  % Magnetometer reference direction (in inertial frame)
std_mag  = simParameters.sensors.mag.std;  % Standard deviation
bias_mag = simParameters.sensors.mag.bias; % Bias of magnetometer
%%% Gyroscope sensor
bias_gyro = simParameters.sensors.gyro.bias;  % Gyroscope bias
std_gyro  = simParameters.sensors.gyro.std;   % Standard deviation

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
    % Generate random star vectors in the inertial frame
    theta = 2*pi*rand(1, number_of_stars);
    phi = acos(2*rand(1, number_of_stars) - 1);
    x_I = sin(phi) .* cos(theta); y_I = sin(phi) .* sin(theta); z_I = cos(phi);
    stars_I = [x_I; y_I; z_I];
    bias_star = simParameters.sensors.star.bias;
    std_star  = simParameters.sensors.star.std;
else
    number_of_stars = 0;
    stars_I = [];
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
    R_k = blkdiag(diag(std_acc_ekf.^2), diag(std_mag_ekf.^2));
end
%%% Actuator parameters
number_of_rw = simParameters.rw.number;
W  = simParameters.rw.W;
%% 2. Simulation containers
x = NaN(7,n); u = NaN(3,n); dq = NaN(4,n-1); T_winf_nosat = NaN(number_of_rw,n);
g_B = NaN(3,n-1); m_B = NaN(3,n-1); stars_B = NaN(3*number_of_stars,n-1);
omega_meas = NaN(3,n-1);
x_est = NaN(7,n-1); y_est = NaN(6+3*number_of_stars,n-3); 
o = NaN(1,n-1); 
%% 3. Initial conditions
% Set initial true state, control, and estimated state
x(:,1)=[simParameters.initialValues.q0; simParameters.initialValues.Wo];
u(:,1) = zeros(3,1);  x_est(:,1) = [1, zeros(1,6)]'; 
T_winf_nosat(:,1) = zeros(number_of_rw, 1);

% Initialize containers for the last known measurement/calculation (Zero-Order Hold).
last_omega_meas = simParameters.initialValues.Wo;
last_u = zeros(3,1);
last_T_winf_nosat = zeros(number_of_rw, 1);

% Get the initial rotation matrix from the true state.
R_init = my_quat2rot(x(1:4, 1));

% Simulate an initial, sensor reading. 
g_B(:,1) = get_accelerometer_reading(R_init, g_I, std_acc, bias_acc); 
m_B(:,1) = get_magnetometer_reading(R_init, m_I, std_mag, bias_mag); 
if simParameters.sensors.star.enable == 1
    stars_B(:,1) = get_star_tracker_reading(R_init, stars_I, std_star, bias_star, number_of_stars);
end

%% 4. Previous calculations to optimize the simulation
Wd_dot = diff(wd')'./diff(t);
hWaitbar = waitbar(0, 'Progress: 0%','Name', 'Attitude simulation progress..');
%% 5. Initialize simulation
for i = 1:n-1
    %% 5.1. Sensor models (add noise, bias and SAMPLING)
    % Gyroscope model with sampling.
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
        if simParameters.ekf.enable == 1
            % Use the pre-calculated, noise-free initial measurements for the TRIAD algorithm.
            q_est_init = triad_algorithm(g_I, m_I, g_B(:,1), m_B(:,1));
            x_est(1:4, i) = q_est_init;
        else
            x_est(1:4, i) = x(1:4,i);
        end
    end
    
    % --- EKF Prediction Step (always runs at high frequency) ---
    % The prediction step runs at each time step 'dt' using the
    % most recent gyroscope measurement (which may be a held measurement).
    A = [eye(4), -0.5 *dt* xi_matrix(x_est(1:4, i)); zeros(3, 4), eye(3)];
    B = 1/2 * dt * [xi_matrix(x_est(1:4, i)); zeros(3, 3)];
    x_pred = A * x_est(:, i) + B * omega_meas(:, i);
    
    Q = dt*[0.25*xi_matrix(x_est(1:4,i))*Q_gyro*xi_matrix(x_est(1:4,i))', zeros(4,3);
              zeros(3,4), Q_gyro];
    P_pred = A * P_cov_ant * A' + Q; 

    % --- EKF Correction Step (runs only when attitude measurements are available) ---
    if mod(i-1, steps_attitude) == 0
        % Generate new attitude sensor measurements using helper functions
        R = my_quat2rot(x(1:4, i));
        g_B(:,i) = get_accelerometer_reading(R, g_I, std_acc, bias_acc);
        m_B(:,i) = get_magnetometer_reading(R, m_I, std_mag, bias_mag);
        
        if simParameters.sensors.star.enable == 1 
            stars_B(:,i) = get_star_tracker_reading(R, stars_I, std_star, bias_star, number_of_stars);
            y = [g_B(:,i); m_B(:,i); stars_B(:,i)];
        else
            y = [g_B(:,i); m_B(:,i)];
        end
        
        % Perform the correction
        C = measurement_jacobian(x_pred(1:4), g_I, m_I, stars_I);
        y_est_pred = C * [x_pred(1:4); zeros(3,1)]; % Predicted measurement (only depends on quaternion)
        y_est(:,i) = y_est_pred;
        k_K = P_pred * C' / (C * P_pred * C' + R_k);  % Kalman gain
        x_est(:, i + 1) = x_pred + k_K * (y - y_est(:,i)); % Correct state
        P_cov = (eye(7) - k_K * C) * P_pred; % Correct covariance
    else
        % If no new measurement, the next state is simply the predicted state
        x_est(:, i + 1) = x_pred;
        P_cov = P_pred;
    end
    
    % Normalize quaternion and update covariance for the next iteration
    x_est(1:4, i + 1) = x_est(1:4, i + 1) / norm(x_est(1:4, i + 1));
    if x_est(1, i + 1) < 0, x_est(1:4, i + 1) = -x_est(1:4, i + 1); end
    P_cov_ant = P_cov;
   
    
    %% 5.3. Control Law (runs at its own sampling rate Ts_control)
    if mod(i-1, steps_control) == 0
        % --- Calculate new control command ---
        if simParameters.ekf.enable == 1
            %%% Use sensors model as feedback signal
            feed_est = [x_est(1:4,i); omega_meas(:,i)];
        else
            %%% Use ideal signals as feedback
            feed_est = x(:, i);
        end
        
        dq(:, i) = Error_quaternio(qd(:, i), feed_est(1:4));
        
        tic
        if simParameters.controller.selector == 1
            u_new = ControlFeedback_rw(I, feed_est, dq(:, i), wd(:, i), Wd_dot(:, i), P, K); 
        else
            k_dot = Gain_estimator_bosk(feed_est(5:7), wd(:, i), dq(:, i), delta, gamma, k, Umax);
            k = k_ant + dt / 6 * (k_dot_ant + 2 * (k_dot_ant + k_dot) + k_dot);
            u_new = Boskovic_control(feed_est(5:7), wd(:, i), dq(:, i), delta, k, Umax);
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
    else
        % --- Hold previous control command (ZOH) ---
        o(i) = NaN; % No control computation cost
    end
    
    % Log the control signal for every step (it will be piecewise constant)
    u(:, i + 1) = last_u;
    T_winf_nosat(:, i + 1) = last_T_winf_nosat;

    %% 5.4. Runge-Kutta 4th order integration for state update
    T_u = W*last_T_winf_nosat;
    g1 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i));
    g2 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g1);
    g3 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g2);
    g4 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g3);
    x(:, i + 1) = x(:, i) + (1 / 6) * (g1 + 2 * g2 + 2 * g3 + g4);
    
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
sensors.meas = [omega_meas; fillmissing([g_B; m_B; stars_B], 'previous',2)];
sensors.est = fillmissing(y_est, 'previous',2);

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
        error_vec = x(1:4,:) - x_est(1:4,:);
        indicators.RMSE = sqrt(mean(error_vec.^2,2));
        %indicators.RSME = cumtrapz(t,error_vec.^2');
         
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

function g_B_meas = get_accelerometer_reading(R, g_I, std_acc, bias_acc)
% Simulates an accelerometer reading.
    W_acc = std_acc .* randn(3, 1);
    g_B_meas = R' * g_I + W_acc + bias_acc;
    g_B_meas = g_B_meas / norm(g_B_meas);
end

function m_B_meas = get_magnetometer_reading(R, m_I, std_mag, bias_mag)
% Simulates a magnetometer reading.
    W_mag = std_mag .* randn(3, 1);
    m_B_meas = R' * m_I + W_mag + bias_mag;
    m_B_meas = m_B_meas / norm(m_B_meas);
end

function stars_B_meas = get_star_tracker_reading(R, stars_I, std_star, bias_star, num_stars)
% Simulates a star tracker reading for multiple stars.
    W_star = std_star .* randn(3, num_stars);
    stars_B_temp = R' * stars_I + W_star + repmat(bias_star, 1, num_stars);
    stars_B_temp = stars_B_temp ./ vecnorm(stars_B_temp, 2);
    stars_B_meas = stars_B_temp(:); % Reshape into a column vector
end