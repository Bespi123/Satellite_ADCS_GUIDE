import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.forces.ForceModel;
import org.orekit.time.AbsoluteDate;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import java.util.ArrayList;
import java.util.List;

// NO MÁS IMPORTACIONES DEL PROPAGADOR AQUÍ

public class StepStorageHandler implements OrekitFixedStepHandler {

    // --- NUEVOS CAMPOS ---
    private final List<ForceModel> forceModels; // Almacenará la lista de fuerzas
    private final AbsoluteDate initialDate;
    private final ArrayList<double[]> history = new ArrayList<>();
    
    // --- NUEVO CONSTRUCTOR ---
    // Ahora recibe la lista de fuerzas directamente desde MATLAB.
    public StepStorageHandler(final List<ForceModel> forceModels, final AbsoluteDate initialDate) {
        this.forceModels = forceModels;
        this.initialDate = initialDate;
    }

    @Override
    public void handleStep(SpacecraftState currentState) {
        double time = currentState.getDate().durationFrom(this.initialDate);
        Vector3D pos = currentState.getPosition();

        double[] stepData = new double[19]; 
        stepData[0] = time;
        stepData[1] = pos.getX();
        stepData[2] = pos.getY();
        stepData[3] = pos.getZ();
        
        // --- USA LA LISTA ALMACENADA DIRECTAMENTE ---
        // Ya no hay necesidad de llamar a propagator.getForceModels()
        for (int i = 0; i < this.forceModels.size(); i++) {
            ForceModel force = this.forceModels.get(i);
            Vector3D acceleration = force.acceleration(currentState, force.getParameters());
            
            int baseIndex = 4 + (i * 3);
            if (baseIndex + 2 < stepData.length) {
                stepData[baseIndex]     = acceleration.getX();
                stepData[baseIndex + 1] = acceleration.getY();
                stepData[baseIndex + 2] = acceleration.getZ();
            }
        }
        
        history.add(stepData);
    }

    public double[][] getHistory() {
        double[][] historyArray = new double[history.size()][19];
        for (int i = 0; i < history.size(); i++) {
            historyArray[i] = history.get(i);
        }
        return historyArray;
    }
}