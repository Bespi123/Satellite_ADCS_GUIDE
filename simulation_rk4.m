function [x, u, T_winf_nosat, x_est, indicators, sensors, error_flag] = simulation_rk4(app,disturbances,simParameters,time)
% SIMULATION_RK4 Simulates CubeSat attitude using a 4th-order Runge-Kutta integrator and EKF.
%
% This program simulates the attitude dynamics of a CubeSat, incorporating
% a 4th-order Runge-Kutta (RK4) integrator for state propagation, an
% Extended Kalman Filter (EKF) for state estimation, and models for
% various sensors and control laws.
%
% Author: bespi123
% Creation Date: [Creation Date, e.g., 2023-10-27]
% Last Modified: [Last Modification Date, e.g., 2024-07-02]
%
% Inputs:
%   app           - UI application object (contains GUI properties like selected controllers).
%   disturbances  - Vector of external disturbance torques applied to the CubeSat over time.
%                   Format: [3xN], where N is the number of time steps.
%   simParameters - Structure containing all simulation parameters, including:
%                   - .setPoint: Desired set points for attitude and angular velocity.
%                   - .initialValues: CubeSat's initial conditions (quaternion, angular velocity, inertia).
%                   - .feedback: Feedback controller gains.
%                   - .boskController: Boskovic controller parameters.
%                   - .sensors: Sensor model parameters (standard deviations, biases).
%                   - .ekf: EKF-specific parameters (standard deviations for the filter model).
%   time          - Structure defining simulation time parameters:
%                   - .t: Time vector [1xN].
%                   - .n: Total number of time steps.
%
% Outputs:
%   x             - Actual CubeSat state over time (quaternion and angular velocity).
%                   Format: [7xN], where rows 1-4 are the quaternion and 5-7 are the angular velocity.
%   u             - Control torque applied by the actuator over time.
%                   Format: [3xN].
%   x_est         - EKF estimated CubeSat state over time (quaternion and gyroscope biases).
%                   Format: [7xN], where rows 1-4 are the estimated quaternion and 5-7 are the biases.
%   indicators    - Structure containing simulation performance metrics and indicators:
%                   - .eulerInt: Cumulative integral of the Euler error angle over time.
%                   - .ascct: Accumulated actuator energy consumption over time.
%                   - .ts: Settlement time for attitude.
%                   - .o: Control law computation time at each step.
%   sensors       - Structure containing sensor data:
%                   - .meas: Simulated sensor measurements (gyroscope, accelerometer, magnetometer, stars).
%                   - .est: EKF predicted measurements.
%   error_flag    - Error flag: 0 if simulation completed without errors, 1 if NaN was detected.
%
% Detailed Description:
%   The `simulation_rk4` function orchestrates the CubeSat attitude dynamics
%   simulation. It operates in a discrete time loop, where at each step:
%   1. Simulated sensor readings with noise and biases are generated.
%   2. An Extended Kalman Filter (EKF) is executed to estimate the current
%      state (orientation and gyroscope biases) from the noisy measurements.
%      The EKF consists of a prediction phase and an update (correction) phase.
%   3. A control law (feedback or Boskovic) is applied based on the attitude
%      error to calculate the necessary control torque.
%   4. The CubeSat's state equation is integrated using the 4th-order
%      Runge-Kutta (RK4) method to propagate the true CubeSat state to the next
%      time step, considering both control and disturbance torques.
%   5. Performance metrics such as quaternion error and actuator consumption are calculated.
%
%   The simulation incorporates models for:
%   - **Accelerometer:** Measures the direction of gravity in the body frame.
%   - **Magnetometer:** Measures the direction of the magnetic field in the body frame.
%   - **Gyroscope:** Measures the angular velocity of the body.
%   - **Star Sensor (optional):** Measures the direction of known stars in the body frame.
%
%   The EKF estimation is based on a 7-element state model (quaternion and
%   gyroscope biases).
%
% Dependencies:
%   - `cubeSatEquationState.m`: Function defining the CubeSat attitude dynamics.
%   - `Error_quaternio.m`: Function to calculate the quaternion error.
%   - `ControlFeedback_rw.m`: Implementation of the feedback controller.
%   - `Boskovic_control.m`: Implementation of the Boskovic controller.
%   - `Gain_estimator_bosk.m`: Gain estimator for the Boskovic controller.
%   - `xi_matrix.m`: Auxiliary function for quaternion operations.
%   - `my_quat2rot.m`: Converts a quaternion to a rotation matrix.
%   - `measurement_jacobian.m`: Calculates the measurement Jacobian matrix for the EKF.
%   - `quat2axang.m` (MATLAB built-in): Converts quaternion to axis-angle.
%   - `quat2eul.m` (MATLAB built-in): Converts quaternion to Euler angles.
%   - `vecnorm.m` (MATLAB built-in): Calculates the norm of vectors.
%   - `cumtrapz.m` (MATLAB built-in): Cumulative numerical integration (trapezoidal method).
%   - `calculateSettlementTime.m`: Function to calculate the settlement time.
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
    %%% Create a diagonal covariance matrix for the star sensor noise
    star = diag(std_star_ekf.^2);  % Diagonal matrix with squared standard deviations
    
    %%% Use the Kronecker product to create a block diagonal matrix for the star sensor
    star_large = kron(eye(number_of_stars), star);  % Block diagonal matrix for multiple stars
    
    %%% Combine accelerometer, magnetometer, and star sensor noise covariance into a single matrix
    R_k = blkdiag(diag(std_acc_ekf.^2), diag(std_mag_ekf.^2), star_large);
else
    %%% If star sensor is not enabled, just combine accelerometer and magnetometer noise covariance
    R_k = blkdiag(diag(std_acc_ekf.^2), diag(std_mag_ekf.^2));
end

%%% Actuator parameters
% Piramidal configuration 1
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

%% 4. Previous calculations to optimize the simulation
%%%Calculate Wd_dot
Wd_dot = diff(wd')'./diff(t);

%%%Create wait bar
hWaitbar = waitbar(0, 'Progress: 0%','Name', 'Attitude simulation progress..');

%% 5. Initialize simulation
for i = 1:n-1
    %% 5.1. Sensor models (add noise and bias)
    R = my_quat2rot(x(1:4, i));  % Rotation matrix from quaternion
    
    % Accelerometer readings (add noise and bias)
    W_acc = std_acc .* randn(3, 1);           % Random noise
    g_B(:, i) = R' * g_I + W_acc + bias_acc;  % Accelerometer measurement
    g_B(:, i) = g_B(:, i)/norm(g_B(:, i));    % Get unitary vector

    % Magnetometer readings (add noise and bias)
    W_mag = std_mag .* randn(3, 1);           % Random noise
    m_B(:, i) = R' * m_I + W_mag + bias_mag;  % Magnetometer measurement
    m_B(:, i) = m_B(:, i)/norm(m_B(:, i));    % Get unitary vector

    % Gyroscope readings (add noise and bias)
    W_gyro = std_gyro .* randn(3, 1);                  % Random noise
    omega_meas(:, i) = x(5:7, i) + W_gyro + bias_gyro; % Gyroscope measurement
    
    %%% Star sensor model
    if simParameters.sensors.star.enable == 1 
        W_star = std_star .* rand(3,number_of_stars);  % Random noise

        % Calculate the star measurements in the body frame (stars_B)
        stars_B_temp = R' * stars_I + W_star + bias_star;  
        % Normalize the measured star vectors
        stars_B_temp = stars_B_temp ./ vecnorm(stars_B_temp, 2);
    
        % Reshape the current star measurements (stars_B_now) into a column vector
        stars_B(:,i) = stars_B_temp(:);  % Reshape into a column vector for further processing or use
        
        % Measurement vector with stars (including accelerometer, magnetometer, and star measurements)
        y = [g_B(:, i); m_B(:, i); stars_B(:,i)];  % Combine accelerometer, magnetometer, and star data into a single measurement vector
    
    else
        % Measurement vector without stars (only accelerometer and magnetometer data)
        y = [g_B(:, i); m_B(:, i)];  % Combine accelerometer and magnetometer data into the measurement vector
    end

    %% 5.2. EKF - Extended Kalman Filter
    if i == 1
        % Initial value calculated using the triad algorithm
        q_est_init = triad_algorithm(g_I, m_I, g_B(:, i), m_B(:, i));
        x_est(1:4, i) = q_est_init;
    end
    % Update the state estimation and control input
    A = [eye(4), -0.5 *dt* xi_matrix(x_est(1:4, i)); zeros(3, 4), eye(3)];
    B = 1/2 * dt * [xi_matrix(x_est(1:4, i)); zeros(3, 3)];
    x_est(:, i + 1) = A * x_est(:, i) + B * omega_meas(:, i);
    
    % Process noise covariance update
    Q = dt*[0.25*xi_matrix(x_est(1:4,i))*Q_gyro*xi_matrix(x_est(1:4,i))', zeros(4,3);
              zeros(3,4), Q_gyro];
    P_cov = A * P_cov_ant * A' + Q; 
    
    % Measurement update (correction step)
    C = measurement_jacobian(x_est(1:4, i), g_I, m_I, stars_I);
    y_est(:,i) = C * x_est(:, i + 1);  % Predicted measurements
    y_est(:,i) = y_est(:,i)/norm(y_est(:,i));

    k_K = P_cov * C' / (C * P_cov * C' + R_k);  % Kalman gain

    % Update state estimate
    x_est(:, i + 1) = x_est(:, i + 1) + k_K * (y - y_est(:,i));

    % Normalize quaternion
    x_est(1:4, i + 1) = x_est(1:4, i + 1) / norm(x_est(1:4, i + 1));
    
    if x_est(1, i + 1) < 0
         x_est(1:4, i + 1) = -x_est(1:4, i + 1);
    end

    % Update covariance
    P_cov = (eye(7) - k_K * C) * P_cov;
    P_cov_ant = P_cov;
    
    %% 5.3. Runge-Kutta 4th order integration for state update
    %T_u = u(:,i);
    T_u = W*T_winf_nosat(:, i);

    g1 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i));
    g2 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g1);
    g3 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g2);
    g4 = dt * cubeSatEquationState(Td(:, i), I, T_u, x(:, i) + 0.5 * g3);
    x(:, i + 1) = x(:, i) + (1 / 6) * (g1 + 2 * g2 + 2 * g3 + g4);
    
    % Calculate quaternion error
    dq(:, i) = Error_quaternio(qd(:, i), x(1:4, i));
    
    % Control law and computational cost estimation
    tic
    if app.controller_popupmenu.ValueIndex == 1
        u(:, i + 1) = ControlFeedback_rw(I, x(:, i), dq(:, i), wd(:, i), Wd_dot(:, i), P, K); 
    else
        k_dot = Gain_estimator_bosk(x(5:7, i), wd(:, i), dq(:, i), delta, gamma, k, Umax);
        % Second-order Simpson integration for k
        k = k_ant + dt / 6 * (k_dot_ant + 2 * (k_dot_ant + k_dot) + k_dot);
        u(:, i + 1) = Boskovic_control(x(5:7, i), wd(:, i), dq(:, i), delta, k, Umax);
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

    figure()
    plot(t,T_winf_nosat);
end

%% 6. Program Functions

function Xi = xi_matrix(q)
% xi_matrix calculates a 3x3 matrix from a 4-element vector q.
%
% Input:
%   q (1x4 vector) - A 4-element vector (typically a quaternion)
%
% Output:
%   Xi (3x3 matrix) - The 3x3 matrix computed from the components of q

% Create the matrix Xi based on the components of q
% Developed by bespi123
    Xi = [ -q(2), -q(3), -q(4);
            q(1), -q(4),  q(3);
            q(4),  q(1), -q(2);
           -q(3),  q(2),  q(1) ];
end

function R = my_quat2rot(q)
% my_quat2rot Converts a quaternion to a 3x3 rotation matrix.
%
%   R = my_quat2rot(q) takes a quaternion `q = [q0; q1; q2; q3]` (scalar
%   component first) and returns the corresponding 3x3 rotation matrix `R`.
%
%   Input:
%       q - A 4x1 vector representing the quaternion `[q_w; q_x; q_y; q_z]`.
%
%   Output:
%       R - A 3x3 rotation matrix.
%
%   The input quaternion is normalized to ensure `R` is a valid rotation matrix.
%   Developed by bespi123

    % Normalize the quaternion
    % Ensures the quaternion is a unit quaternion for a valid rotation.
    q = q / norm(q);

    % Extract components
    q0 = q(1); % Scalar part
    q1 = q(2); % Vector x-component
    q2 = q(3); % Vector y-component
    q3 = q(4); % Vector z-component

    % Construct rotation matrix
    % Standard formula for quaternion to rotation matrix conversion.
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end

function C = measurement_jacobian(q, g_I, m_I, star_I)
% This function calculates the measurement Jacobian matrix for an Extended Kalman Filter (EKF)
% used in attitude estimation. The Jacobian relates changes in the state
% vector (quaternion and biases) to changes in the measured vectors (gravity,
% magnetic field, and star vectors).
%
% Inputs:
%   q       - Quaternion representing the current attitude [q0, q1, q2, q3] (scalar first).
%   g_I     - Gravity vector in the inertial frame [gx, gy, gz].
%   m_I     - Magnetic field vector in the inertial frame [mx, my, mz].
%   star_I  - Matrix of star vectors in the inertial frame. Each column
%             represents a star vector: [star1_x, star2_x, ...;
%                                       star1_y, star2_y, ...;
%                                       star1_z, star2_z, ...].
%
% Output:
%   C       - 6x7 (or larger, depending on the number of stars) measurement
%             Jacobian matrix. The first four columns correspond to the
%             quaternion components, and the remaining three columns are zeros
%             representing the partial derivatives with respect to biases
%             (which are assumed to not directly affect the measurement
%             equations in this formulation).
% Developed by bespi123

    % Extract quaternion components
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    % Extract gravity components
    gx = g_I(1); gy = g_I(2); gz = g_I(3);
    % Extract magnetic field components
    rx = m_I(1); ry = m_I(2); rz = m_I(3);
    % Create the Jacobian matrix for gravity measurement
    C_g = 2*[gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, -gx*q2 + gy*q1 - gz*q0, -gx*q3 + gy*q0 + gz*q1, 0, 0, 0;
            -gx*q3 + gy*q0 + gz*q1,  gx*q2 - gy*q1 + gz*q0,  gx*q1 + gy*q2 + gz*q3, -gx*q0 - gy*q3 + gz*q2, 0, 0, 0;
             gx*q2 - gy*q1 + gz*q0,  gx*q3 - gy*q0 - gz*q1,  gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, 0, 0, 0];
    % Create the Jacobian matrix for magnetic field measurement
    C_m = 2*[rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, -rx*q2 + ry*q1 - rz*q0, -rx*q3 + ry*q0 + rz*q1, 0, 0, 0;
            -rx*q3 + ry*q0 + rz*q1,  rx*q2 - ry*q1 + rz*q0,  rx*q1 + ry*q2 + rz*q3, -rx*q0 - ry*q3 + rz*q2, 0, 0, 0;
             rx*q2 - ry*q1 + rz*q0,  rx*q3 - ry*q0 - rz*q1,  rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, 0, 0, 0];
    
    % Initialize C_star matrix based on the number of star vectors
    stars_num = size(star_I,2);
    C_star = NaN(3*stars_num,7);
    
    % Calculate Jacobian for each star vector
    for i=1:stars_num
        star_x = star_I(1,i);
        star_y = star_I(2,i);
        star_z = star_I(3,i);
        
        % Define row indices for the current star's Jacobian contribution
        index1 = ((i-1)*3+1);
        index2 = (index1+3)-1;
        
        % Calculate the 3x7 Jacobian for the current star
        C_star(index1:index2,:) = 2*[star_x*q0 + star_y*q3 - star_z*q2,  star_x*q1 + star_y*q2 + star_z*q3, -star_x*q2 + star_y*q1 - star_z*q0, -star_x*q3 + star_y*q0 + star_z*q1, 0, 0, 0;
            -star_x*q3 + star_y*q0 + star_z*q1,  star_x*q2 - star_y*q1 + star_z*q0,  star_x*q1 + star_y*q2 + star_z*q3, -star_x*q0 - star_y*q3 + star_z*q2, 0, 0, 0;
             star_x*q2 - star_y*q1 + star_z*q0,  star_x*q3 - star_y*q0 - star_z*q1,  star_x*q0 + star_y*q3 - star_z*q2,  star_x*q1 + star_y*q2 + star_z*q3, 0, 0, 0];
    end
    
    % Concatenate all Jacobian contributions to form the final measurement Jacobian matrix
    C = [C_g; C_m; C_star];
end

function q_triad = triad_algorithm(r1,r2,b1,b2)
% triad_algorithm Calculates the attitude matrix using the TRIAD algorithm.
%
%   A_triad = triad_algorithm(r1, r2, b1, b2) computes the attitude matrix
%   (also known as the rotation matrix) that describes the orientation
%   of a rigid body. It achieves this by utilizing two non-collinear
%   vector observations, one pair from a known reference frame (e.g., inertial)
%   and the corresponding pair observed from the body frame whose attitude
%   is to be determined.
%
%   Inputs:
%   r1 - 3x1 (or compatible) reference vector for the first observation.
%        This vector should be expressed in the reference frame.
%   r2 - 3x1 (or compatible) reference vector for the second observation.
%        This vector must be non-collinear with 'r1'.
%   b1 - 3x1 (or compatible) body vector corresponding to 'r1'.
%        This is the same physical vector as 'r1', but expressed in the
%        body frame.
%   b2 - 3x1 (or compatible) body vector corresponding to 'r2'.
%        This is the same physical vector as 'r2', but expressed in the
%        body frame, and must be non-collinear with 'b1'.
%
%   Output:
%   q_triad - A 4x1 attitude (rotation) quaternion.

    %%% Ensure input vectors are unitary column vectors
    r1 = reshape(r1, [], 1); 
    r2 = reshape(r2, [], 1); 
    b1 = reshape(b1, [], 1); 
    b2 = reshape(b2, [], 1);
    
    %%% Construct an orthonormal basis (triad) in the reference frame (V-frame)
    %v1 = r1;           % Aligned with the first reference vector.
    v2 = cross(r1,r2) / norm(cross(r1,r2)); % Perpendicular to the plane formed by r1 and r2
    %v3 = cross(v1,v2); % Completes the right-handed orthonormal triad 
    
    %%% Construct an orthonormal basis (triad) in the body frame (W-frame)
    %w1 = b1;           % Aligned with the first body vector.
    w2 = cross(b1,b2) / norm(cross(b1,b2)); % Perpendicular to the plane formed by b1 and b2
    %w3 = cross(w1,w2); % Completes the right-handed orthonormal triad
    
    %%% Estimate the Attitude Matrix (A_triad)
    % A more common (and often more intuitive) formulation is:
    % A_triad = [w1, w2, w3] * [v1, v2, v3]';
    % Both forms are mathematically equivalent and should yield the same result.
    A_triad = b1*r1' + (cross(b1,w2))*(cross(r1,v2))' + w2*v2';

    q_triad = dcm2quat(A_triad)';
end