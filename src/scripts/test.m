% 1. CONFIGURACIÓN INICIAL
clear; clc;

% Asegúrate de que MATLAB apunte al entorno de Python correcto
% (donde instalaste 'ahrs'). Cambia la ruta si es necesario.
%pyenv("Version", "D:\git\Satellite_ADCS_GUIDE\.venv\Scripts\python.exe"); 

% --- 2. GENERAR DATOS DE SENSORES SIMULADOS (en MATLAB) ---
num_samples = 1000;
sample_rate = 100.0; % 100 Hz
t = (0:num_samples-1)' / sample_rate; % Vector columna

% Datos del giroscopio (rotación en eje Z)
gyr_data = zeros(num_samples, 3);
gyr_data(:, 3) = sin(2 * pi * 0.1 * t) * pi; % rad/s

% Datos del acelerómetro (estático con ruido)
acc_data = zeros(num_samples, 3);
acc_data(:, 3) = 9.81; % Gravedad
acc_data = acc_data + randn(num_samples, 3) * 0.1; % Ruido

% --- 3. PROCESAMIENTO CON EKF LLAMADO DESDE MATLAB ---
disp('Iniciando el procesamiento llamando a Python desde MATLAB...');

% Crear una instancia del objeto EKF de Python
ekf = py.ahrs.filters.EKF(pyargs('frequency', sample_rate));

% Pre-alocar el array de cuaterniones en MATLAB (usando 1-based indexing)
Q = zeros(num_samples, 4); 


% Estimación inicial (OJO: MATLAB usa paréntesis y empieza en el índice 1)
Q(1, :) = double(py.ahrs.common.orientation.acc2q(acc_data(1, :)));

% Bucle para procesar todas las muestras
for t = 2:num_samples
    % --- ¡AQUÍ ESTÁ LA LLAMADA DIRECTA! ---
    % Pasamos arrays de MATLAB directamente al método .update del objeto de Python.
    Q(t, :) = double(ekf.update(Q(t-1, :), gyr_data(t, :), acc_data(t, :)));
end

disp('Procesamiento completado.');
disp('Forma del array de cuaterniones resultante:');
disp(size(Q));
disp('Ejemplo del último cuaternión calculado:');
disp(Q(end, :));