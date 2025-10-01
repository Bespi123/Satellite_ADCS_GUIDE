classdef StarTracker
    % STARTRACKER Encapsula el modelo de un sensor de estrellas.
    % Esta es una clase por valor, ya que no necesita mantener un estado 
    % interno que cambie con el tiempo (a diferencia del IMU con deriva).

    properties
        enable          logical % Habilita (true) o deshabilita (false) el sensor
        numberOfStars   double  % Número de estrellas a simular
        std             double  % Desviación estándar del ruido (vector 3x1)
        bias            double  % Sesgo del sensor (vector 3x1)
        stars_I         double  % Matriz 3xN con los vectores de dirección de las estrellas en el marco inercial
    end
    
    methods
        function obj = StarTracker(star_params)
            % CONSTRUCTOR: Inicializa las propiedades del sensor de estrellas
            % a partir de una estructura de parámetros.
            
            obj.enable = star_params.enable;
            
            if obj.enable
                obj.numberOfStars = star_params.numberOfStars;
                obj.std = star_params.std;
                obj.bias = star_params.bias;
                
                % Generar y almacenar los vectores de estrellas aleatorios
                % para que sean consistentes durante toda la simulación.
                num_stars = obj.numberOfStars;
                theta = 2 * pi * rand(1, num_stars);
                phi = acos(2 * rand(1, num_stars) - 1);
                x_I = sin(phi) .* cos(theta);
                y_I = sin(phi) .* sin(theta);
                z_I = cos(phi);
                obj.stars_I = [x_I; y_I; z_I];
            else
                % Si está deshabilitado, establece valores vacíos o nulos.
                obj.numberOfStars = 0;
                obj.std = [];
                obj.bias = [];
                obj.stars_I = [];
            end
        end
        
        function reading = getReading(obj, R)
            % GETREADING Simula una medición del sensor de estrellas.
            % 'R' es la matriz de rotación del marco inercial al de cuerpo (I->B).
            
            % Si el sensor está apagado, no devuelve ninguna lectura.
            if ~obj.enable
                reading = [];
                return;
            end
            
            % 1. Transformar los vectores de estrellas al marco del cuerpo.
            stars_B_true = R' * obj.stars_I;
            
            % 2. Generar ruido blanco gaussiano para cada estrella.
            num_stars = obj.numberOfStars;
            white_noise = obj.std .* randn(3, num_stars);
            
            % 3. Añadir el ruido y el bias (replicado para cada estrella).
            bias_rep = repmat(obj.bias, 1, num_stars);
            stars_B_measured = stars_B_true + white_noise + bias_rep;
            
            % 4. Normalizar cada vector, ya que un sensor de estrellas mide direcciones.
            stars_B_normalized = stars_B_measured ./ vecnorm(stars_B_measured, 2);
            
            % 5. Aplanar la matriz de salida a un único vector columna para el EKF.
            reading = stars_B_normalized(:);
        end
    end
end