classdef AttitudeEKF < handle
    % ATTITUDEEKF Encapsula la lógica de un Extended Kalman Filter para estimación de actitud.
    % Mantiene el estado estimado (x_est) y la covarianza (P_cov) internamente.

    properties
        x_est       double % Vector de estado estimado [q; w] (7x1)
        P_cov       double % Matriz de covarianza del error (7x7)
        
        Q_gyro      double % Matriz de covarianza del ruido del proceso (giroscopio) (3x3)
        R           double % Matriz de covarianza del ruido de la medición (NxN)
        
        % Referencias a los objetos de los sensores
        imu         % Objeto de la clase IMU
        starTracker % Objeto de la clase StarTracker
    end
    
    methods
        function obj = AttitudeEKF(ekf_params, imu_obj, st_obj)
            % CONSTRUCTOR: Inicializa el EKF con sus parámetros y los objetos de los sensores.
            
            % Guardar referencias a los objetos de sensores
            obj.imu = imu_obj;
            obj.starTracker = st_obj;
            
            % Inicializar estado y covarianza
            obj.x_est = [1; 0; 0; 0; 0; 0; 0]; % Estado inicial por defecto (sin rotación)
            obj.P_cov = eye(7) * 1e-2; % Incertidumbre inicial
            
            % Configurar matriz de ruido del proceso (Q)
            obj.Q_gyro = diag(ekf_params.gyro.std.^2);
            
            % Configurar matriz de ruido de la medición (R)
            R_acc = diag(ekf_params.acc.std.^2);
            R_mag = diag(ekf_params.mag.std.^2);
            if obj.starTracker.enable
                num_stars = obj.starTracker.numberOfStars;
                R_star_single = diag(ekf_params.star.std.^2);
                R_star_large = kron(eye(num_stars), R_star_single);
                obj.R = blkdiag(R_acc, R_mag, R_star_large);
            else
                obj.R = blkdiag(R_acc, R_mag);
            end
        end
        
        function initializeWithTriad(obj, g_B, m_B)
            % Inicializa la actitud del EKF usando el algoritmo TRIAD.
            q_triad = obj.triad_algorithm(obj.imu.g_I, obj.imu.m_I, g_B, m_B);
            obj.x_est(1:4) = q_triad;
        end
        
        function predict(obj, omega_meas_filtered, dt)
            % Realiza el paso de predicción del EKF.
            q = obj.x_est(1:4);
            
            % Matriz de transición de estado A
            A = [eye(4), -0.5 * dt * obj.xi_matrix(q); 
                 zeros(3, 4), eye(3)];
            
            % Matriz de entrada B
            B = 0.5 * dt * [obj.xi_matrix(q); zeros(3, 3)];
            
            % Predicción del estado
            obj.x_est = A * obj.x_est + B * omega_meas_filtered;
            
            % Matriz de ruido del proceso Q (dependiente del estado)
            Xi_q = obj.xi_matrix(q);
            Q_top_left = 0.25 * Xi_q * obj.Q_gyro * Xi_q';
            Qd = dt * blkdiag(Q_top_left, obj.Q_gyro);
            
            % Predicción de la covarianza
            obj.P_cov = A * obj.P_cov * A' + Qd;
        end
        
        function correct(obj, y_meas)
            % Realiza el paso de corrección del EKF.
            
            % Calcular el Jacobiano de la medición
            C = obj.computeMeasurementJacobian();
            
            % Calcular la ganancia de Kalman
            K_gain = obj.P_cov * C' / (C * obj.P_cov * C' + obj.R);
            
            % Calcular la medición predicha
            y_pred = obj.predictMeasurement();
            
            % Corregir el estado
            obj.x_est = obj.x_est + K_gain * (y_meas - y_pred);
            
            % Corregir la covarianza
            obj.P_cov = (eye(7) - K_gain * C) * obj.P_cov;
            
            % Normalizar el cuaternio y asegurar q0 positivo
            obj.x_est(1:4) = obj.x_est(1:4) / norm(obj.x_est(1:4));
            if obj.x_est(1) < 0
                obj.x_est(1:4) = -obj.x_est(1:4);
            end
        end
    end
    
    methods (Access = private)
        % Métodos de ayuda encapsulados dentro de la clase
        
        function y_pred = predictMeasurement(obj)
            % Predice la medición de los sensores basada en el estado actual.
            R_pred = obj.quat2rot(obj.x_est(1:4));
            g_B_pred = R_pred' * obj.imu.g_I;
            m_B_pred = R_pred' * obj.imu.m_I;
            
            if obj.starTracker.enable
                stars_B_pred_matrix = R_pred' * obj.starTracker.stars_I;
                stars_B_pred = stars_B_pred_matrix(:);
                y_pred = [g_B_pred; m_B_pred; stars_B_pred];
            else
                y_pred = [g_B_pred; m_B_pred];
            end
        end

        function C = computeMeasurementJacobian(obj)
            % Calcula el Jacobiano C.
            q = obj.x_est(1:4);
            q0=q(1); q1=q(2); q2=q(3); q3=q(4);
            
            % Jacobiano para el acelerómetro
            g_I = obj.imu.g_I;
            gx=g_I(1); gy=g_I(2); gz=g_I(3);
            C_g = 2*[gx*q0+gy*q3-gz*q2, gx*q1+gy*q2+gz*q3, -gx*q2+gy*q1-gz*q0, -gx*q3+gy*q0+gz*q1;
                     -gx*q3+gy*q0+gz*q1, gx*q2-gy*q1+gz*q0, gx*q1+gy*q2+gz*q3, -gx*q0-gy*q3+gz*q2;
                     gx*q2-gy*q1+gz*q0, gx*q3-gy*q0-gz*q1, gx*q0+gy*q3-gz*q2, gx*q1+gy*q2+gz*q3];
            
            % Jacobiano para el magnetómetro
            m_I = obj.imu.m_I;
            rx=m_I(1); ry=m_I(2); rz=m_I(3);
            C_m = 2*[rx*q0+ry*q3-rz*q2, rx*q1+ry*q2+rz*q3, -rx*q2+ry*q1-rz*q0, -rx*q3+ry*q0+rz*q1;
                     -rx*q3+ry*q0+rz*q1, rx*q2-ry*q1+rz*q0, rx*q1+ry*q2+rz*q3, -rx*q0-ry*q3+rz*q2;
                     rx*q2-ry*q1+rz*q0, rx*q3-ry*q0-rz*q1, rx*q0+ry*q3-rz*q2, rx*q1+ry*q2+rz*q3];
            
            C_star = [];
            if obj.starTracker.enable
                star_I = obj.starTracker.stars_I;
                num_stars = size(star_I, 2);
                C_star = zeros(3*num_stars, 4);
                for i=1:num_stars
                    sx=star_I(1,i); sy=star_I(2,i); sz=star_I(3,i);
                    C_star_i = 2*[sx*q0+sy*q3-sz*q2, sx*q1+sy*q2+sz*q3, -sx*q2+sy*q1-sz*q0, -sx*q3+sy*q0+sz*q1;
                                 -sx*q3+sy*q0+sz*q1, sx*q2-sy*q1+sz*q0, sx*q1+sy*q2+sz*q3, -sx*q0-sy*q3+sz*q2;
                                 sx*q2-sy*q1+sz*q0, sx*q3-sy*q0-sz*q1, sx*q0+sy*q3-sz*q2, sx*q1+sy*q2+sz*q3];
                    C_star(3*(i-1)+1:3*i, :) = C_star_i;
                end
            end
            
            C_attitude = [C_g; C_m; C_star];
            C = [C_attitude, zeros(size(C_attitude,1), 3)]; % Jacobiano completo (7 columnas)
        end
        
        function R = quat2rot(~, q)
            q = q / norm(q);
            q0=q(1); q1=q(2); q2=q(3); q3=q(4);
            R = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
                 2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
                 2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
        end

        function Xi = xi_matrix(~, q)
            Xi = [-q(2), -q(3), -q(4);
                   q(1), -q(4),  q(3);
                   q(4),  q(1), -q(2);
                  -q(3),  q(2),  q(1)];
        end

        function q_triad = triad_algorithm(~, r1, r2, b1, b2)
            r1=r1(:); r2=r2(:); b1=b1(:); b2=b2(:);
            t1_r = r1;
            t1_b = b1;
            t2_r = cross(r1,r2)/norm(cross(r1,r2));
            t2_b = cross(b1,b2)/norm(cross(b1,b2));
            t3_r = cross(t1_r, t2_r);
            t3_b = cross(t1_b, t2_b);
            A_r = [t1_r, t2_r, t3_r];
            A_b = [t1_b, t2_b, t3_b];
            A_triad = A_b * A_r';
            
            % Conversión de Matriz de Rotación a Cuaternio
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