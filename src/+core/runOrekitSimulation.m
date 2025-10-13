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

end % End of function