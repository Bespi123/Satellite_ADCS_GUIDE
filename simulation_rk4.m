function [x, u, x_est, indicators, sensors, error_flag] = simulation_rk4(app,disturbances,simParameters,time)
%%%Program developed by: bespi123

%%Recover variables
t = time.t; n=time.n;
qd = repmat(simParameters.setPoint.qd',n, 1)';    %%Desired attitude array
wd = repmat(simParameters.setPoint.Wd',n, 1)';    %%Desired angular rate array 
Td = disturbances; I = simParameters.initialValues.I;

%%%Controller gains
if app.controller_popupmenu.ValueIndex == 1
    P = simParameters.feedback.Peye; K = simParameters.feedback.Peye;
else
    delta = simParameters.boskController.delta;
    gamma = simParameters.boskController.gamma;
    k     = simParameters.boskController.k0;
    Umax  = simParameters.boskController.Umax;

    k_ant = k; k_dot_ant = 0;
end

%%% Sensors model
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
P_cov_ant = eye(7);  % Initial state covariance matrix (7x7 identity matrix)

% If star sensor is enabled, include the star sensor noise in the covariance matrix
if simParameters.sensors.star.enable == 1
    % Create a diagonal covariance matrix for the star sensor noise
    star = diag(std_star.^2);  % Diagonal matrix with squared standard deviations
    
    % Use the Kronecker product to create a block diagonal matrix for the star sensor
    star_large = kron(eye(number_of_stars), star);  % Block diagonal matrix for multiple stars
    
    % Combine accelerometer, magnetometer, and star sensor noise covariance into a single matrix
    R_k = blkdiag(diag(std_acc.^2), diag(std_mag.^2), star_large);
else
    % If star sensor is not enabled, just combine accelerometer and magnetometer noise covariance
    R_k = blkdiag(diag(std_acc.^2), diag(std_mag.^2));
end

%%% Error fl
error_flag = 0;

%%%Simulation containers
x = NaN(7,n); u = NaN(3,n); dq = NaN(4,n-1);
g_B = NaN(3,n-1); m_B = NaN(3,n-1); stars_B = NaN(3*number_of_stars,n-1);
omega_meas = NaN(3,n-1);
x_est = NaN(7,n); y_est = NaN(6+3*number_of_stars,n-3); 
o = NaN(1,n-1); 

%%%eulerInt = zeros(1,n); asscct = zeros(1,n); 
%%%Initial conditions
x(:,1)=[simParameters.initialValues.q0; simParameters.initialValues.Wo];
u(:,1) = zeros(3,1);  x_est(:,1) = [1, zeros(1,6)]'; 

%%%Calculate Wd_dot
Wd_dot = diff(wd')'./diff(t);

%%%Create wait bar
%%%Progress bar settings
hWaitbar = waitbar(0, 'Progress: 0%','Name', 'Attitude simulation progress..');

for i = 1:n-1
    % Calculate time step
    dt = t(2) - t(1);  
    
    % Sensor models (add noise and bias)
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


    % EKF - Extended Kalman Filter

    % Update the state estimation and control input
    A = [eye(4), -0.5 *dt* xi_matrix(x_est(1:4, i)); zeros(3, 4), eye(3)];
    B = 1/2 * dt * [xi_matrix(x_est(1:4, i)); zeros(3, 3)];
    
    x_est(:, i + 1) = A * x_est(:, i) + B * omega_meas(:, i);
    
    % Process noise covariance update
    Q_gyro = diag(std_gyro.^2);
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
    
    % Runge-Kutta 4th order integration for state update
    g1 = dt * cubeSatEquationState(Td(:, i), I, u(:, i), x(:, i));
    g2 = dt * cubeSatEquationState(Td(:, i), I, u(:, i), x(:, i) + 0.5 * g1);
    g3 = dt * cubeSatEquationState(Td(:, i), I, u(:, i), x(:, i) + 0.5 * g2);
    g4 = dt * cubeSatEquationState(Td(:, i), I, u(:, i), x(:, i) + 0.5 * g3);
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

    %%
    if error_flag == 0
        %%%EULERINT calculation
        ang_and_axis = quat2axang(dq'); 
        eulerang = ang_and_axis(:,4);
        indicators.eulerInt = cumtrapz(t(1:end-1),eulerang); 
        %%%ASCCT Calculation
        ascct_dot=vecnorm(u,2,1).^2;
        indicators.ascct = cumtrapz(t,ascct_dot); 
        %%%Settlement time calculation
        tol = 5/100; % 2%
        indicators.ts = calculateSettlementTime(180/pi*quat2eul(dq'), t, tol);
        %%%Computation time
        indicators.o = o;
    end
end

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
           q(1), -q(4), q(3);
           q(4), q(1), -q(2);
          -q(3), q(2), q(1) ];
end



function R = my_quat2rot(q)
    % Normalizar el cuaternión primero
    q = q / norm(q);

    % Extraer componentes
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    % Construir matriz de rotación
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end

function C = measurement_jacobian(q, g_I, m_I, star_I)
% Measurement Jacobian matrix for EKF
% Inputs:
%   q - quaternion [q0, q1, q2, q3] (scalar first)
%   g_I - gravity vector in inertial frame [gx, gy, gz]
%   m_I - magnetic field vector in inertial frame [mx, my, mz]
%
% Output:
%   C - 6x7 measurement Jacobian matrix

% Extract quaternion components
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

% Extract gravity components
gx = g_I(1); gy = g_I(2); gz = g_I(3);

% Extract magnetic field components
rx = m_I(1); ry = m_I(2); rz = m_I(3);

% Create the Jacobian matrix
C_g = 2*[gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, -gx*q2 + gy*q1 - gz*q0, -gx*q3 + gy*q0 + gz*q1, 0, 0, 0;
        -gx*q3 + gy*q0 + gz*q1,  gx*q2 - gy*q1 + gz*q0,  gx*q1 + gy*q2 + gz*q3, -gx*q0 - gy*q3 + gz*q2, 0, 0, 0;
         gx*q2 - gy*q1 + gz*q0,  gx*q3 - gy*q0 - gz*q1,  gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, 0, 0, 0];
C_m = 2*[rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, -rx*q2 + ry*q1 - rz*q0, -rx*q3 + ry*q0 + rz*q1, 0, 0, 0;
        -rx*q3 + ry*q0 + rz*q1,  rx*q2 - ry*q1 + rz*q0,  rx*q1 + ry*q2 + rz*q3, -rx*q0 - ry*q3 + rz*q2, 0, 0, 0;
         rx*q2 - ry*q1 + rz*q0,  rx*q3 - ry*q0 - rz*q1,  rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, 0, 0, 0];

stars_num = size(star_I,2);
C_star = NaN(3*stars_num,7);

for i=1:stars_num
    star_x = star_I(1,i);
    star_y = star_I(2,i);
    star_z = star_I(3,i);
    
    index1 = ((i-1)*3+1);
    index2 = (index1+3)-1;

    C_star(index1:index2,:) = 2*[star_x*q0 + star_y*q3 - star_z*q2,  star_x*q1 + star_y*q2 + star_z*q3, -star_x*q2 + star_y*q1 - star_z*q0, -star_x*q3 + star_y*q0 + star_z*q1, 0, 0, 0;
        -star_x*q3 + star_y*q0 + star_z*q1,  star_x*q2 - star_y*q1 + star_z*q0,  star_x*q1 + star_y*q2 + star_z*q3, -star_x*q0 - star_y*q3 + star_z*q2, 0, 0, 0;
         star_x*q2 - star_y*q1 + star_z*q0,  star_x*q3 - star_y*q0 - star_z*q1,  star_x*q0 + star_y*q3 - star_z*q2,  star_x*q1 + star_y*q2 + star_z*q3, 0, 0, 0];
end

C = [C_g; C_m; C_star];
end