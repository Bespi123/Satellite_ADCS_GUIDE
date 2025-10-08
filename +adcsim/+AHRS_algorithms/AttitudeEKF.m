classdef AttitudeEKF < handle
    % ATTITUDEEKF Encapsulates the logic of an Extended Kalman Filter for attitude estimation.
    % It maintains the estimated state (x_est) and covariance (P_cov) internally.
    properties
        x_est       double % Estimated state vector [q; gyro_bias] (7x1)
        P_cov       double % Error covariance matrix (7x7)
        
        Q_gyro      double % Process noise covariance matrix (gyroscope) (3x3)
        R           double % Measurement noise covariance matrix (NxN)
        
        dt          double % EKF sample time

        % References to sensor objects
        imu         % IMU class object
        starTracker % StarTracker class object
    end
    
    methods (Access = public)
        function obj = AttitudeEKF(ekf_params, sample_time, imu_obj, st_obj)
            import adcsim.utils.*
                        
            % Save references to sensor objects
            obj.imu = imu_obj;
            obj.dt  = sample_time;
            % Check if a StarTracker was provided
            if nargin > 3 && ~isempty(st_obj)
                obj.starTracker = st_obj;
            else
                obj.starTracker = struct('enable', false); % Default to disabled
            end
            
            % Initialize state and covariance
            obj.x_est = [1; 0; 0; 0; 0; 0; 0]; % Default initial state (no rotation)
            obj.P_cov = eye(7) * 1e-2; % Initial uncertainty
            
            % Configure process noise matrix (Q)
            % Q_gyro represents the noise for the 3 error state elements corresponding to angular velocity.
            obj.Q_gyro = diag(ekf_params.gyro.std.^2);
            
            % Configure measurement noise matrix (R)
            R_acc = diag(ekf_params.acc.std.^2);
            R_mag = diag(ekf_params.mag.std.^2);
            
            % Combine R matrices based on available sensors
            if obj.starTracker.enable
                % Note: R_star assumes star measurement noise parameters are included in ekf_params
                num_stars = obj.starTracker.numberOfStars;
                R_star_single = diag(ekf_params.star.std.^2);
                R_star_large = kron(eye(num_stars), R_star_single); % Stack R for N stars
                obj.R = blkdiag(R_acc, R_mag, R_star_large); % Combined R matrix: [R_acc, R_mag, R_star]
            else
                obj.R = blkdiag(R_acc, R_mag); % Combined R matrix: [R_acc, R_mag]
            end
        end
        
        function initializeWithTriad(obj, g_B, m_B)
            % Initializes the EKF attitude using the TRIAD algorithm.
            q_triad = obj.triad_algorithm(obj.imu.g_I, obj.imu.m_I, g_B, m_B);
            obj.x_est(1:4) = q_triad;
        end
        
        function predict(obj, gyroscope)
            % Performs the EKF prediction step.
            q = obj.x_est(1:4);
            
            % State transition matrix A (linearized for small time step)
            % A = [d(q)/d(q), d(q)/d(w); d(w)/d(q), d(w)/d(w)]
            A = [eye(4), -0.5 * obj.dt * obj.xi_matrix(q); 
                 zeros(3, 4), eye(3)];
            
            % Input matrix B (for the measured angular velocity input)
            B = 0.5 * obj.dt * [obj.xi_matrix(q); zeros(3, 3)];
            
            % State prediction: x_est(k+1) = A*x_est(k) + B*u(k)
            % Note: This is an *approximate* prediction, common in discrete EKF
            obj.x_est = A * obj.x_est + B * gyroscope;
            
            % Process noise covariance matrix Qd (dependent on state)
            % Qd = dt * G * Q_gyro * G' where G = [0.5*Xi_q; I] in the error state formulation
            Xi_q = obj.xi_matrix(q);
            Q_top_left = 0.25 * Xi_q * obj.Q_gyro * Xi_q';
            Qd = obj.dt * blkdiag(Q_top_left, obj.Q_gyro);
            
            % Covariance prediction: P_cov(k+1) = A*P_cov(k)*A' + Qd
            obj.P_cov = A * obj.P_cov * A' + Qd;
        end
        
        function correct(obj, y_meas)
            % Performs the EKF correction step.
            
            % Calculate the Measurement Jacobian (H or C)
            C = obj.computeMeasurementJacobian();
            
            % Calculate the Kalman Gain: K = P*C' * (C*P*C' + R)^(-1)
            K_gain = obj.P_cov * C' / (C * obj.P_cov * C' + obj.R);
            
            % Calculate the predicted measurement (h(x))
            y_pred = obj.predictMeasurement();
            
            % Correct the state: x_est = x_est + K*(y_meas - y_pred)
            obj.x_est = obj.x_est + K_gain * (y_meas - y_pred);
            
            % Correct the covariance: P_cov = (I - K*C) * P_cov
            obj.P_cov = (eye(7) - K_gain * C) * obj.P_cov;
            
            % Normalize the quaternion and ensure positive q0 (for convention)
            obj.x_est(1:4) = obj.x_est(1:4) / norm(obj.x_est(1:4));
            if obj.x_est(1) < 0
                obj.x_est(1:4) = -obj.x_est(1:4);
            end
        end
    end
    
    methods (Access = private)
        % Helper methods encapsulated within the class
        
        function y_pred = predictMeasurement(obj)
            % Predicts the sensor measurement based on the current state.
            q_est = obj.x_est(1:4);
            R_pred = obj.quat2rot(q_est); % Rotation matrix from Body to Inertial (C_I^B)
            
            % Predicted measurements for IMU references in the Body frame
            g_B_pred = R_pred' * obj.imu.g_I; % R' is C_B^I
            m_B_pred = R_pred' * obj.imu.m_I;
            
            if obj.starTracker.enable
                % Predicted StarTracker measurements (stars in Body frame)
                stars_B_pred_matrix = R_pred' * obj.starTracker.stars_I;
                stars_B_pred = stars_B_pred_matrix(:); % Flatten into a vector
                y_pred = [g_B_pred; m_B_pred; stars_B_pred];
            else
                y_pred = [g_B_pred; m_B_pred];
            end
        end
        
        function C = computeMeasurementJacobian(obj)
            % Calculates the Measurement Jacobian C (or H).
            q = obj.x_est(1:4);
            q0=q(1); q1=q(2); q2=q(3); q3=q(4);
            
            % Jacobian for the accelerometer (derivative of h_g = R' * g_I w.r.t q)
            g_I = obj.imu.g_I;
            gx=g_I(1); gy=g_I(2); gz=g_I(3);
            % C_g is a 3x4 matrix
            C_g = 2*[gx*q0+gy*q3-gz*q2, gx*q1+gy*q2+gz*q3, -gx*q2+gy*q1-gz*q0, -gx*q3+gy*q0+gz*q1;
                     -gx*q3+gy*q0+gz*q1, gx*q2-gy*q1+gz*q0, gx*q1+gy*q2+gz*q3, -gx*q0-gy*q3+gz*q2;
                     gx*q2-gy*q1+gz*q0, gx*q3-gy*q0-gz*q1, gx*q0+gy*q3-gz*q2, gx*q1+gy*q2+gz*q3];
            
            % Jacobian for the magnetometer (derivative of h_m = R' * m_I w.r.t q)
            m_I = obj.imu.m_I;
            rx=m_I(1); ry=m_I(2); rz=m_I(3);
            % C_m is a 3x4 matrix
            C_m = 2*[rx*q0+ry*q3-rz*q2, rx*q1+ry*q2+rz*q3, -rx*q2+ry*q1-rz*q0, -rx*q3+ry*q0+rz*q1;
                     -rx*q3+ry*q0+rz*q1, rx*q2-ry*q1+rz*q0, rx*q1+ry*q2+rz*q3, -rx*q0-ry*q3+rz*q2;
                     rx*q2-ry*q1+rz*q0, rx*q3-ry*q0-rz*q1, rx*q0+ry*q3-rz*q2, rx*q1+ry*q2+rz*q3];
            
            C_star = [];
            if obj.starTracker.enable
                % Jacobian for the StarTracker
                star_I = obj.starTracker.stars_I;
                num_stars = size(star_I, 2);
                C_star = zeros(3*num_stars, 4);
                for i=1:num_stars
                    % Jacobian C_star_i is a 3x4 matrix for the i-th star
                    sx=star_I(1,i); sy=star_I(2,i); sz=star_I(3,i);
                    C_star_i = 2*[sx*q0+sy*q3-sz*q2, sx*q1+sy*q2+sz*q3, -sx*q2+sy*q1-sz*q0, -sx*q3+sy*q0+sz*q1;
                                 -sx*q3+sy*q0+sz*q1, sx*q2-sy*q1+sz*q0, sx*q1+sy*q2+sz*q3, -sx*q0-sy*q3+sz*q2;
                                 sx*q2-sy*q1+sz*q0, sx*q3-sy*q0-sz*q1, sx*q0+sy*q3-sz*q2, sx*q1+sy*q2+sz*q3];
                    C_star(3*(i-1)+1:3*i, :) = C_star_i; % Stack vertically
                end
            end
            
            % Combine attitude Jacobians: C_attitude is [C_g; C_m; C_star]
            C_attitude = [C_g; C_m; C_star]; 
            
            % Full Jacobian C: [C_attitude, zeros(N_meas, 3)] (N_meas x 7 matrix)
            % The last 3 columns are for the angular velocity error (d/dw), which is zero.
            C = [C_attitude, zeros(size(C_attitude,1), 3)]; 
        end
        
        function R_mat = quat2rot(~, q)
            % Converts a quaternion to a 3x3 rotation matrix (C_I^B).
            q = q / norm(q);
            q0=q(1); q1=q(2); q2=q(3); q3=q(4);
            R_mat = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
                 2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
                 2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
        end
        
        function Xi = xi_matrix(~, q)
            % Computes the Xi matrix used in quaternion derivative/propagation.
            % Xi(q) * w = q_dot, where q_dot is the derivative w.r.t. the angular velocity vector.
            Xi = [-q(2), -q(3), -q(4);
                   q(1), -q(4),  q(3);
                   q(4),  q(1), -q(2);
                  -q(3),  q(2),  q(1)];
        end
        
        function q_triad = triad_algorithm(~, r1, r2, b1, b2)
            % Implements the TRIAD algorithm for initial attitude estimation.
            % r1, r2: Reference vectors in Inertial frame.
            % b1, b2: Measured vectors in Body frame.
            r1=r1(:); r2=r2(:); b1=b1(:); b2=b2(:);
            
            % Create the orthogonal triads
            t1_r = r1;
            t1_b = b1;
            t2_r = cross(r1,r2)/norm(cross(r1,r2));
            t2_b = cross(b1,b2)/norm(cross(b1,b2));
            t3_r = cross(t1_r, t2_r);
            t3_b = cross(t1_b, t2_b);
            
            % Assemble the matrices and compute the rotation matrix C_B^I (R_mat = A_b * A_r')
            A_r = [t1_r, t2_r, t3_r];
            A_b = [t1_b, t2_b, t3_b];
            A_triad = A_b * A_r';
            
            % Convert Rotation Matrix to Quaternion (standard algorithm)
            tr = trace(A_triad);
            if tr > 0
                S = sqrt(tr+1.0) * 2;
                qw = 0.25 * S;
                qx = (A_triad(3,2) - A_triad(2,3)) / S;
                qy = (A_triad(1,3) - A_triad(3,1)) / S;
                qz = (A_triad(2,1) - A_triad(1,2)) / S;
            elseif (A_triad(1,1) > A_triad(2,2)) && (A_triad(1,1) > A_triad(3,3))
                S = sqrt(1.0 + A_triad(1,1) - A_triad(2,2) - A_triad(3,3)) * 2;
                qw = (A_triad(3,2) - A_triad(2,3)) / S;
                qx = 0.25 * S;
                qy = (A_triad(1,2) + A_triad(2,1)) / S;
                qz = (A_triad(1,3) + A_triad(3,1)) / S;
            elseif A_triad(2,2) > A_triad(3,3)
                S = sqrt(1.0 + A_triad(2,2) - A_triad(1,1) - A_triad(3,3)) * 2;
                qw = (A_triad(1,3) - A_triad(3,1)) / S;
                qx = (A_triad(1,2) + A_triad(2,1)) / S;
                qy = 0.25 * S;
                qz = (A_triad(2,3) + A_triad(3,2)) / S;
            else
                S = sqrt(1.0 + A_triad(3,3) - A_triad(1,1) - A_triad(2,2)) * 2;
                qw = (A_triad(2,1) - A_triad(1,2)) / S;
                qx = (A_triad(1,3) + A_triad(3,1)) / S;
                qy = (A_triad(2,3) + A_triad(3,2)) / S;
                qz = 0.25 * S;
            end
            q_triad = [qw; qx; qy; qz];
        end
    end
end