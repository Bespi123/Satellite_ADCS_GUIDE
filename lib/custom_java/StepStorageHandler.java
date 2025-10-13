import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.forces.ForceModel;
import org.orekit.time.AbsoluteDate;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import java.util.ArrayList;
import java.util.List;

/**
 * A custom Orekit step handler designed to store detailed simulation data for analysis in MATLAB.
 * <p>
 * This class is called by an Orekit propagator at fixed time intervals. At each step, it records:
 * <ul>
 * <li>The simulation time elapsed since the initial date.</li>
 * <li>The spacecraft's position vector (x, y, z).</li>
 * <li>The acceleration vector produced by <strong>each individual force model</strong>.</li>
 * </ul>
 * The data is stored in an ArrayList and can be retrieved as a 2D double array, which is ideal for MATLAB.
 */

public class StepStorageHandler implements OrekitFixedStepHandler {

    // --- CLASS FIELDS ---
    /** The list of force models used in the simulation, provided by the calling environment (MATLAB). */
    private final List<ForceModel> forceModels;
    /** The start date of the propagation, used as the reference for time stamps (t=0). */
    private final AbsoluteDate initialDate;
    /** The main storage for the simulation data. Each entry is a double array representing one time step. */
    private final ArrayList<double[]> history = new ArrayList<>();
    
    /**
     * Constructs the step handler.
     * <p>
     * This constructor receives all necessary dependencies from the calling environment (MATLAB)
     * before the propagation starts.
     *
     * @param forceModels The list of {@link ForceModel} objects whose accelerations will be recorded.
     * @param initialDate The absolute start date of the simulation.
     */
    public StepStorageHandler(final List<ForceModel> forceModels, final AbsoluteDate initialDate) {
        this.forceModels = forceModels;
        this.initialDate = initialDate;
    }

    /**
     * This method is called by the Orekit propagator at each fixed time step.
     * <p>
     * It extracts the current state, calculates the acceleration from each individual force model,
     * and stores all the data in the history list.
     *
     * @param currentState The current state of the spacecraft provided by the propagator.
     */
    @Override
    public void handleStep(SpacecraftState currentState) {
        // Get the elapsed time in seconds from the start of the simulation.
        double time = currentState.getDate().durationFrom(this.initialDate);
        // Get the current position vector
        Vector3D pos = currentState.getPosition();
        
        // Prepare a new array to hold all data for this specific time step.
        // Array size is 19: 1 (time) + 3 (position) + 5 forces * 3 (acceleration components) = 19
        double[] stepData = new double[19]; 
        stepData[0] = time;
        stepData[1] = pos.getX();
        stepData[2] = pos.getY();
        stepData[3] = pos.getZ();
        
        // --- Use the stored list of forces directly ---
        // There is no need to call propagator.getForceModels() anymore.
        for (int i = 0; i < this.forceModels.size(); i++) {
            ForceModel force = this.forceModels.get(i);
            // Calculate the acceleration vector for this specific force model.
            Vector3D acceleration = force.acceleration(currentState, force.getParameters());
            
            // Calculate the base index in the array for this force's data.
            // Force 0 (Gravity) -> index 4; Force 1 (Drag) -> index 7, etc.
            int baseIndex = 4 + (i * 3);
            if (baseIndex + 2 < stepData.length) {
                stepData[baseIndex]     = acceleration.getX();
                stepData[baseIndex + 1] = acceleration.getY();
                stepData[baseIndex + 2] = acceleration.getZ();
            }
        }

        // Add the complete data array for this step to the history list.
        history.add(stepData);
    }

    /**
     * Retrieves the complete simulation history as a 2D array.
     * <p>
     * This method is intended to be called from MATLAB after the propagation is complete.
     * It converts the internal ArrayList into a primitive 2D double array, which is
     * easily handled by MATLAB.
     *
     * @return A 2D double array where each row represents a time step and the columns
     * contain the stored data [time, x, y, z, ax1, ay1, az1, ax2, ...].
     */
    public double[][] getHistory() {
        double[][] historyArray = new double[history.size()][19];
        for (int i = 0; i < history.size(); i++) {
            historyArray[i] = history.get(i);
        }
        return historyArray;
    }
}