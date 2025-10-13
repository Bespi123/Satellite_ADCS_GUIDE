function history = runOrekitSimulation(simParams)
    %RUNOREKITSIMULATION Runs a high-fidelity orbit propagation using Orekit.
    %
    %   history = RUNOREKITSIMULATION(simParams) propagates the orbit of a
    %   spacecraft based on the parameters defined in the simParams struct.
    %
    %   INPUTS:
    %       simParams - A MATLAB struct with the following fields:
    %         .initialDate   - MATLAB datetime object for the simulation start.
    %         .initialOrbit  - Struct with .a, .e, .i, .raan, .pa, .ta (Keplerian elements).
    %         .spacecraft    - Struct with .mass, .dragArea, .dragCoeff, .srpArea, .srpCoeff.
    %         .duration      - Struct with .num_orbits.
    %         .integrator    - Struct with .minStep, .maxStep, .posTolerance.
    %         .output        - Struct with .stepSize (in seconds for data logging).
    %
    %   OUTPUTS:
    %       history - A matrix containing the time history of the simulation. Each
    %                 row corresponds to a time step and has 19 columns:
    %                 [time, pos_x, pos_y, pos_z, acc_grav_x, ..., acc_srp_z]
    %
    %   Requires the StepStorageHandler.class to be compiled and available.
    
    %% 2. Import Orekit Java Classes
    % --- Propagation & Integration ---
    import org.orekit.propagation.numerical.NumericalPropagator;
    import org.orekit.propagation.SpacecraftState;
    import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
    % --- Time & Coordinates ---
    import org.orekit.frames.FramesFactory;
    import org.orekit.time.AbsoluteDate;
    import org.orekit.time.TimeScalesFactory;
    import org.orekit.orbits.PositionAngleType;
    import org.orekit.orbits.KeplerianOrbit;
    % --- Force Models & Bodies ---
    import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
    import org.orekit.forces.gravity.ThirdBodyAttraction;
    import org.orekit.forces.radiation.SolarRadiationPressure;
    import org.orekit.forces.radiation.IsotropicRadiationSingleCoefficient;
    import org.orekit.forces.drag.DragForce;
    import org.orekit.forces.drag.IsotropicDrag;
    import org.orekit.utils.Constants;
    import org.orekit.bodies.CelestialBodyFactory;
    import java.util.ArrayList;
    
    %% STEP 1: Configure Numerical Integrator
    disp('Configuring numerical propagator...');
    integrator = DormandPrince853Integrator(simParams.integrator.minStep, ...
                                            simParams.integrator.maxStep, ...
                                            simParams.integrator.posTolerance, ...
                                            simParams.integrator.posTolerance);
    
    %% STEP 2: Define Initial State
    utc = TimeScalesFactory.getUTC();
    % Convert MATLAB datetime to Orekit AbsoluteDate
    dt = simParams.initialDate;
    initialDate = AbsoluteDate(dt.Year, dt.Month, dt.Day, dt.Hour, dt.Minute, dt.Second, utc);
    
    inertialFrame = FramesFactory.getGCRF();
    p = simParams.initialOrbit; % Shortcut for parameters
    initialOrbit = KeplerianOrbit(p.a, p.e, p.i, p.pa, p.raan, p.ta, ...
                                  PositionAngleType.TRUE, inertialFrame, initialDate, Constants.WGS84_EARTH_MU);
    initialState = SpacecraftState(initialOrbit, simParams.spacecraft.mass);
    
    %% STEP 3: Create Propagator
    propagator = NumericalPropagator(integrator);
    
    %% STEP 4: Add Detailed Force Models
    force_models_list = ArrayList();
    sc = simParams.spacecraft; % Shortcut for spacecraft parameters
    
    % --- Force 1: Earth Gravity (21x21 spherical harmonics) ---
    gravityProvider = org.orekit.forces.gravity.potential.GravityFieldFactory.getNormalizedProvider(21, 21);
    gravityForce = HolmesFeatherstoneAttractionModel(inertialFrame, gravityProvider);
    propagator.addForceModel(gravityForce);
    force_models_list.add(gravityForce);
    
    % --- Force 2: Atmospheric Drag ---
    sun = CelestialBodyFactory.getSun();
    earth = org.orekit.bodies.OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, ...
        Constants.WGS84_EARTH_FLATTENING, FramesFactory.getITRF(org.orekit.utils.IERSConventions.IERS_2010, true));
    atmosphere = org.orekit.models.earth.atmosphere.HarrisPriester(sun, earth);
    dragForce = DragForce(atmosphere, IsotropicDrag(sc.dragArea, sc.dragCoeff));
    propagator.addForceModel(dragForce);
    force_models_list.add(dragForce);
    
    % --- Force 3: Third Body Gravity (Sun & Moon) ---
    moon = CelestialBodyFactory.getMoon();
    sunAttraction = ThirdBodyAttraction(sun);
    moonAttraction = ThirdBodyAttraction(moon);
    propagator.addForceModel(sunAttraction);
    force_models_list.add(sunAttraction);
    propagator.addForceModel(moonAttraction);
    force_models_list.add(moonAttraction);
    
    % --- Force 4: Solar Radiation Pressure (SRP) ---
    isotropicRadiation = IsotropicRadiationSingleCoefficient(sc.srpArea, sc.srpCoeff);
    srpForce = SolarRadiationPressure(sun, earth, isotropicRadiation);
    propagator.addForceModel(srpForce);
    force_models_list.add(srpForce);
    
    disp('Propagator created with high-fidelity force models.');
    fprintf('\n');
    
    %% STEP 5: Configure Handler and Run Propagation
    fixedStepHandler = StepStorageHandler(force_models_list, initialDate);
    propagator.setStepHandler(simParams.output.stepSize, fixedStepHandler);
    
    propagation_duration = initialOrbit.getKeplerianPeriod() * simParams.duration.num_orbits;
    disp(['Propagating orbit for ', num2str(simParams.duration.num_orbits), ' orbits...']);
    
    propagator.setInitialState(initialState);
    finalState = propagator.propagate(initialDate.shiftedBy(propagation_duration));
    disp('Propagation complete.');
    fprintf('\n');
    
    %% STEP 6: Retrieve and Return Results
    disp('Retrieving simulation history...');
    history_java = fixedStepHandler.getHistory();
    
    % Convert Java 2D array to a MATLAB matrix
    [rows, cols] = size(history_java);
    history = zeros(rows, cols);
    for r = 1:rows
        for c = 1:cols
            history(r,c) = history_java(r,c);
        end
    end
    disp('Data successfully retrieved.');

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


end % End of function