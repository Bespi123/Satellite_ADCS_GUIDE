classdef IMU < handle
    % IMU Encapsula los sensores de una Unidad de Medición Inercial.
    % Incluye el giroscopio, acelerómetro y magnetómetro.

    properties
        Gyroscope       struct
        Accelerometer   struct
        Magnetometer    struct
        g_I             double % Vector de referencia de gravedad
        m_I             double % Vector de referencia magnético
    end
    
    methods
        function obj = IMU(sensors, g_I, m_I)            
            % Create the structures directly before assigning them

            obj.Gyroscope = struct(...
                'std', sensors.gyro.std, ...
                'bias', sensors.gyro.bias, ...
                'tau_gyro_filter', sensors.gyro.tau, ...
                'last_filtered_omega', zeros(3,1)...
            );
            
            obj.Accelerometer = struct(...
                'std', sensors.acc.std, ...
                'bias', sensors.acc.bias ...
            );
            
            obj.Magnetometer = struct(...
                'std', sensors.mag.std, ...
                'bias', sensors.mag.bias ...
            );

            obj.g_I = g_I;
            obj.m_I = m_I;
        end
        
        % ... (methods remain the same) ...
        
        function reading = getGyroscopeReading(obj, omega)
            noise = obj.Gyroscope.std .* randn(3, 1);
            reading = omega + noise + obj.Gyroscope.bias;
        end

        function reading_filtered = filterGyroscopeReading(obj, gyro_reading, Ts_gyro)       
            % First-Order Low-Pass Filter for Gyroscope Measurement
            alpha = Ts_gyro / (obj.Gyroscope.tau_gyro_filter + Ts_gyro);
            reading_filtered = (1 - alpha) * obj.Gyroscope.last_filtered_omega + alpha * gyro_reading;
            obj.Gyroscope.last_filtered_omega = reading_filtered;
        end
        
        function reading = getAccelerometerReading(obj, R)
            noise = obj.Accelerometer.std .* randn(3, 1);
            reading = R' * obj.g_I + noise + obj.Accelerometer.bias;
            reading = reading / norm(reading);
        end
        
        function reading = getMagnetometerReading(obj, R)
            noise = obj.Magnetometer.std .* randn(3, 1);
            reading = R' * obj.m_I + noise + obj.Magnetometer.bias;
            reading = reading / norm(reading);
        end
    end
end