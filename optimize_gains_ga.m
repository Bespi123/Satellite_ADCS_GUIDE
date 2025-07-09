function optimize_gains_ga(app,disturbances,simParameters,time)
% --- 1. Configuración Inicial ---
% Necesitas una instancia de tu app y los parámetros iniciales.
% Una forma es ejecutar tu app, y luego ejecutar este script desde la consola de MATLAB.
% Asegúrate de que la variable 'app' exista en el workspace.
% if ~exist('app', 'var') || ~isvalid(app)
%     error('Por favor, ejecuta la app ADCS_program primero y asegúrate de que la variable "app" esté en el workspace.');
% end

% % Carga los parámetros iniciales
% simParameters = app.simulationParameters;
% time = struct('ti', app.initTime_edit.Value, 'tf', app.finalTime_edit.Value, 'step', app.step_edit.Value);
% time.t = time.ti:time.step:time.tf;
% time.n = length(time.t);
% disturbances = getDisturbances(app.simulationParameters.disturbances, time); % Asumiendo que getDisturbances está accesible

% --- 2. Configuración del Algoritmo Genético ---
% Define la función objetivo
fitnessFcn = @(gains) objectiveFunction(gains, app, disturbances, simParameters, time);

% Número de variables a optimizar (6 para el controlador Feedback: 3 para P, 3 para K)
nvars = 6; 

% Límites inferiores (Lower Bounds) y superiores (Upper Bounds) para las ganancias
% ¡Estos límites son cruciales! Debes definirlos con sensatez.
lb = [0, 0, 0, 0, 0, 0]; % Valores mínimos para cada ganancia
ub = [10,   10,   10,   10,   10,   10];   % Valores máximos para cada ganancia

% %%% --- DEFINIR LA POBLACIÓN INICIAL --- %%%
% Define aquí uno o más conjuntos de ganancias que creas que son buenos.
% Cada fila es un individuo.
initial_gains_1 = [0.2, 0.2, 0.2, 0.3, 0.3, 0.3]; % Un buen punto de partida conocido
% Crea la matriz que se pasará al GA.
initial_population = [initial_gains_1];

% Opciones del GA
options = optimoptions('ga', ...
    'PopulationSize', 5, ...       % Individuos por generación
    'MaxGenerations', 10, ...      % Número de generaciones
    'Display', 'iter', ...          % Muestra el progreso en cada iteración
    'PlotFcn', {@gaplotbestf, @gaplotstopping}, ... % Gráfica el mejor fitness y el criterio de parada
    'UseParallel', false, ...
    'InitialPopulationMatrix', initial_population);           % ¡MUY RECOMENDADO para acelerar el proceso!

% --- 3. Ejecutar el Algoritmo Genético ---
disp('Iniciando optimización con Algoritmo Genético...');
[optimal_gains, min_cost] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], options);

% --- 4. Mostrar Resultados ---
disp('Optimización finalizada.');
fprintf('Costo mínimo encontrado: %.4f\n', min_cost);
disp('Ganancias óptimas encontradas (P y K diagonales):');
disp(optimal_gains);

% Actualiza los parámetros en tu app con los resultados óptimos
app.simulationParameters.feedback.Peye = diag(optimal_gains(1:3));
app.simulationParameters.feedback.Keye = diag(optimal_gains(4:6));
msgbox('¡Ganancias óptimas encontradas y cargadas en la app!', 'Optimización Completa');