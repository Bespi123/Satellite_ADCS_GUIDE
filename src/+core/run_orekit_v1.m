% numerical_propagation_matlab_SIMULACION.m
% -------------------------------------------------------------------------
% OBJETIVO: Ejecutar la simulación de la órbita con Orekit.
% PRE-REQUISITO: El archivo 'StepStorageHandler.class' ya debe existir
%                en la misma carpeta (compilado previamente).
% -------------------------------------------------------------------------
clear; clc; clear java;

% % %% 1. Configuración del Entorno (SIN COMPILACIÓN)
% % disp('--- SCRIPT DE SIMULACIÓN ---');
% % disp('Configurando el entorno de Orekit...');
% % 
% % % Definir las rutas a las carpetas
% % lib_folder = 'orekit';
% % data_folder = 'orekit/orekit-data';
% % 
% % % --- Añadir todas las librerías .jar al path de Java ---
% % javaaddpath(fullfile(lib_folder, 'orekit-13.0.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-clustering-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-core-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-fft-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-filtering-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-fitting-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-geometry-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-ode-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-optim-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-samples-4.0.2.jar'));
% % javaaddpath(fullfile(lib_folder, 'hipparchus-stat-4.0.2.jar'));
% % 
% % % --- AÑADIR LA CARPETA ACTUAL AL PATH DE JAVA ---
% % % !! ESTE PASO ES CRUCIAL !!
% % % Le dice a MATLAB dónde encontrar el 'StepStorageHandler.class' que ya compilamos.
% % javaaddpath(fullfile(pwd,'myJavaCodes')); 
% % disp('Manejador de pasos (.class) añadido al path.');
% % 
% % % --- Configurar el proveedor de datos de Orekit ---
% % import org.orekit.data.*;
% % import java.io.File;
% % manager = DataContext.getDefault().getDataProvidersManager();
% % manager.clearProviders();
% % manager.addProvider(DirectoryCrawler(File(data_folder)));
% % 
% % disp('Entorno listo.');
% % fprintf('\n');

%% 2. Importar las clases de Orekit necesarias
% --- Integración y Propagación ---
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.SpacecraftState;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
% --- Coordenadas y Tiempo ---
import org.orekit.frames.FramesFactory;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.orbits.PositionAngleType;
import org.orekit.orbits.KeplerianOrbit;
% --- Modelos de Fuerza ---
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.ThirdBodyAttraction;
import org.orekit.forces.radiation.SolarRadiationPressure;
import org.orekit.forces.radiation.IsotropicRadiationSingleCoefficient;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.models.earth.atmosphere.*;
import org.orekit.utils.Constants;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.propagation.sampling.OrekitStepNormalizer;

%% PASO 1: Configurar el Integrador Numérico
disp('Configurando el propagador numérico...');
minStep = 0.001; % segundos
maxStep = 1000;  % segundos
posTolerance = 1.0; % metros
integrator = DormandPrince853Integrator(minStep, maxStep, posTolerance, posTolerance);

%% PASO 2: Definir el Estado Inicial Preciso
utc = TimeScalesFactory.getUTC();
initialDate = AbsoluteDate(2025, 10, 10, 23, 30, 00.000, utc);
inertialFrame = FramesFactory.getGCRF();
a = Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 420e3; % Radio Tierra + 420 km altitud
e = 0.0001;
i = 51.6 * pi / 180; % Inclinación en radianes
raan = 0.0;
pa = 0.0;
ta = 0.0;
mu = Constants.WGS84_EARTH_MU;
initialOrbit = KeplerianOrbit(a, e, i, pa, raan, ta, PositionAngleType.TRUE, inertialFrame, initialDate, mu);
initialMass = 10.0; % kg
initialState = SpacecraftState(initialOrbit, initialMass);

%% PASO 3: Crear el Propagador (inicialmente vacío)
propagator = NumericalPropagator(integrator);

%% PASO 4: Añadir los Modelos de Fuerza Detallados
% --- NUEVO: Crear una lista de Java para guardar nuestras fuerzas ---
import java.util.ArrayList;
force_models_list = ArrayList();

% --- Fuerza 1: Gravedad de la Tierra (con armónicos esféricos 21x21) ---
gravityProvider = org.orekit.forces.gravity.potential.GravityFieldFactory.getNormalizedProvider(21, 21);
gravityForce = HolmesFeatherstoneAttractionModel(inertialFrame, gravityProvider);
propagator.addForceModel(gravityForce);
force_models_list.add(gravityForce); % <-- AÑADIR A NUESTRA LISTA

% --- Fuerza 2: Arrastre Atmosférico (Drag) ---
sun = CelestialBodyFactory.getSun();
earth = org.orekit.bodies.OneAxisEllipsoid( ...
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS, ...
    Constants.WGS84_EARTH_FLATTENING, ...
    FramesFactory.getITRF(org.orekit.utils.IERSConventions.IERS_2010, true));
atmosphere = org.orekit.models.earth.atmosphere.HarrisPriester(sun, earth);
dragCrossSection = 1.0; % m^2
dragCoeff = 2.2;
dragForce = DragForce(atmosphere, IsotropicDrag(dragCrossSection, dragCoeff));
propagator.addForceModel(dragForce);
force_models_list.add(dragForce); % <-- AÑADIR A NUESTRA LISTA

% --- Fuerza 3: Gravedad de Terceros Cuerpos (Sol y Luna) ---
moon = CelestialBodyFactory.getMoon();
sunAttraction = ThirdBodyAttraction(sun);
moonAttraction = ThirdBodyAttraction(moon);
propagator.addForceModel(sunAttraction);
force_models_list.add(sunAttraction); % <-- AÑADIR A NUESTRA LISTA
propagator.addForceModel(moonAttraction);
force_models_list.add(moonAttraction); % <-- AÑADIR A NUESTRA LISTA

% --- Fuerza 4: Presión de Radiación Solar (SRP) ---
srpCrossSection = 1.0; % m^2
cr = 1.2;
isotropicRadiation = org.orekit.forces.radiation.IsotropicRadiationSingleCoefficient(srpCrossSection, cr);
srpForce = org.orekit.forces.radiation.SolarRadiationPressure(sun, earth, isotropicRadiation);
propagator.addForceModel(srpForce);
force_models_list.add(srpForce); % <-- AÑADIR A NUESTRA LISTA

disp('Propagador numérico creado con los siguientes modelos de fuerza:');
disp('- Gravedad terrestre (21x21)');
disp('- Arrastre atmosférico (HarrisPriester)');
disp('- Gravedad del Sol y la Luna');
disp('- Presión de Radiación Solar');
fprintf('\n');

%% PASO 5: Configurar Manejador y Ejecutar la Propagación
% --- CORREGIDO: Pasar nuestra lista autoconstruida al constructor ---
fixedStepHandler = StepStorageHandler(force_models_list, initialDate);

% Establecer el manejador (esta línea sigue igual)
propagator.setStepHandler(60.0, fixedStepHandler);

% Ejecutar la propagación
num_orbits = 3;
propagation_duration = initialOrbit.getKeplerianPeriod() * num_orbits;
disp(['Propagando la órbita por ', num2str(num_orbits), ' órbitas...']);
propagator.setInitialState(initialState);
finalState = propagator.propagate(initialDate.shiftedBy(propagation_duration));
disp('Propagación completada.');
fprintf('\n');

%% 6. Graficar los Resultados
disp('Graficando la trayectoria y las fuerzas individuales...');

% --- OBTENER LOS DATOS ---
% Llama al método correcto para obtener la matriz de historial (19 columnas)
history_java = fixedStepHandler.getHistory();

% Convertir el arreglo 2D de Java a una matriz de MATLAB
[rows, cols] = size(history_java);
history = zeros(rows, cols);
for r = 1:rows
    for c = 1:cols
        history(r,c) = history_java(r,c);
    end
end

% --- GRÁFICO 1: TRAYECTORIA 3D ---
figure; % Crea la primera figura para la órbita
hold on;

% CORREGIDO: Usa las columnas correctas para la posición (ahora son 2, 3 y 4)
plot3(history(:,2)/1000, history(:,3)/1000, history(:,4)/1000, 'b-', 'LineWidth', 2);

% El resto del código para dibujar la Tierra y los marcadores es el mismo
[x, y, z] = sphere(50);
earth_radius_km = Constants.WGS84_EARTH_EQUATORIAL_RADIUS / 1000;
surf(x*earth_radius_km, y*earth_radius_km, z*earth_radius_km, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot3(history(1,2)/1000, history(1,3)/1000, history(1,4)/1000, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot3(history(end,2)/1000, history(end,3)/1000, history(end,4)/1000, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

% Etiquetas y configuración del gráfico
title('Órbita Propagada con Orekit');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
grid on; axis equal;
legend('Trayectoria', 'Tierra', 'Inicio', 'Fin', 'Location', 'best');
view(45, 25);
hold off;


% GRÁFICO 2: FUERZAS INDIVIDUALES

figure; % Crea una segunda figura para las fuerzas
hold on;

% --- Extraer los datos de la matriz ---
% El tiempo está en la columna 1. Lo convertimos a horas para que se vea mejor.
time_vector = history(:, 1) / 3600; 

% Extraer los vectores de aceleración de sus respectivas columnas
accel_gravity = history(:, 5:7);   % Columnas para Gravedad (21x21)
accel_drag    = history(:, 8:10);  % Columnas para Arrastre
accel_sun     = history(:, 11:13); % Columnas para Sol (3er cuerpo)
accel_moon    = history(:, 14:16); % Columnas para Luna (3er cuerpo)
accel_srp     = history(:, 17:19); % Columnas para Presión de Radiación Solar

% --- Calcular la magnitud de cada fuerza ---
% Usamos vecnorm para calcular la magnitud de cada vector fila
mag_gravity = vecnorm(accel_gravity, 2, 2);
mag_drag    = vecnorm(accel_drag, 2, 2);
mag_sun     = vecnorm(accel_sun, 2, 2);
mag_moon    = vecnorm(accel_moon, 2, 2);
mag_srp     = vecnorm(accel_srp, 2, 2);

% --- Graficar en Escala Logarítmica ---
% La gravedad es millones de veces más fuerte que las otras fuerzas.
% Usar semilogy() es la única forma de poder ver y comparar todas en un mismo gráfico.
semilogy(time_vector, mag_gravity, 'LineWidth', 2);
semilogy(time_vector, mag_drag, 'LineWidth', 2);
semilogy(time_vector, mag_sun, 'LineWidth', 2);
semilogy(time_vector, mag_moon, 'LineWidth', 2);
semilogy(time_vector, mag_srp, 'LineWidth', 2);

% Etiquetas y configuración del gráfico
title('Magnitud de las Aceleraciones por Fuerza (Escala Logarítmica)');
xlabel('Tiempo de Simulación (horas)');
ylabel('Aceleración (m/s^2)');
legend('Gravedad (21x21)', 'Arrastre Atmosférico', 'Tercer Cuerpo (Sol)', 'Tercer Cuerpo (Luna)', 'Presión Radiación Solar');
grid on;
hold off;