classdef StarTracker
    %STARTRACKER Encapsulates the model of a star tracker sensor.
    %   This is a value class, as it does not need to maintain an internal
    %   state that changes over time (unlike the IMU with its bias drift).
    %   Date          Author          Notes
    %   08/20/2025    Bespi123        Create simple Star tracker model
    %% Public Properties
    properties
        enable          logical % Enables (true) or disables (false) the sensor.
        numberOfStars   double  % Number of stars to be simulated.
        std             double  % Standard deviation of the measurement noise (3x1 vector).
        bias            double  % Sensor bias (3x1 vector).
        stars_I         double  % 3xN matrix containing the direction vectors of the stars in the inertial frame.
    end
    
    methods
        function obj = StarTracker(star_params)
            % CONSTRUCTOR: Initializes the star tracker's properties
            % from a given parameters structure.
            
            obj.enable = star_params.enable;
            
            if obj.enable
                obj.numberOfStars = star_params.numberOfStars;
                obj.std = star_params.std;
                obj.bias = star_params.bias;
                
                % Generate and store random star vectors so they remain
                % consistent throughout the entire simulation.
                num_stars = obj.numberOfStars;
                theta = 2 * pi * rand(1, num_stars);
                phi = acos(2 * rand(1, num_stars) - 1); % Uniformly distributed points on a sphere
                x_I = sin(phi) .* cos(theta);
                y_I = sin(phi) .* sin(theta);
                z_I = cos(phi);
                obj.stars_I = [x_I; y_I; z_I];
            else
                % If disabled, set properties to empty or null values.
                obj.numberOfStars = 0;
                obj.std = [];
                obj.bias = [];
                obj.stars_I = [];
            end
        end
        
        function reading = getReading(obj, R)
            % GETREADING Simulates a measurement from the star tracker.
            % 'R' is the rotation matrix from the body frame to the inertial frame (B->I).
            
            % If the sensor is disabled, return no reading.
            if ~obj.enable
                reading = [];
                return;
            end
            
            % 1. Transform the star vectors from the inertial frame to the body frame.
            %    We use the transpose of the (B->I) rotation matrix, which is (I->B).
            stars_B_true = R' * obj.stars_I;
            
            % 2. Generate Gaussian white noise for each star.
            num_stars = obj.numberOfStars;
            white_noise = obj.std .* randn(3, num_stars);
            
            % 3. Add the noise and bias (replicated for each star).
            bias_rep = repmat(obj.bias, 1, num_stars);
            stars_B_measured = stars_B_true + white_noise + bias_rep;
            
            % 4. Normalize each vector, as a star tracker measures directions (unit vectors).
            stars_B_normalized = stars_B_measured ./ vecnorm(stars_B_measured, 2);
            
            % 5. Flatten the output matrix into a single column vector for use in the filter (e.g., EKF).
            reading = stars_B_normalized(:);
        end
    end
end