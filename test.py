import numpy as np
from ahrs.filters import EKF
from ahrs.common.orientation import acc2q

# --- 1. Generar datos de sensores simulados (ESTA ES LA PARTE QUE FALTABA) ---
# En un caso real, cargarías estos datos desde un archivo o un sensor.
num_samples = 1000
sample_rate = 100.0  # 100 Hz
t = np.arange(num_samples) / sample_rate

# Datos del giroscopio (simulando una rotación suave en el eje Z)
gyr_data = np.zeros((num_samples, 3))
gyr_data[:, 2] = np.sin(2 * np.pi * 0.1 * t) * np.pi  # rad/s

# Datos del acelerómetro (simulando que está mayormente estático)
# Se añade un poco de ruido para que sea más realista.
acc_data = np.zeros((num_samples, 3))
acc_data[:, 2] = 9.81  # Gravedad en el eje Z
acc_data += np.random.randn(num_samples, 3) * 0.1  # Ruido gaussiano

# --- 2. Tu código para el procesamiento con el EKF (SIN CAMBIOS) ---
print("Iniciando el procesamiento con el Filtro Extendido de Kalman...")

ekf = EKF(frequency=sample_rate) # Es buena práctica indicar la frecuencia
Q = np.zeros((num_samples, 4))
print(  "Forma del array de cuaterniones inicial:", Q.shape)
# Estimación inicial a partir de la primera muestra del acelerómetro
Q[0] = acc2q(acc_data[0])
print("Cuaternión inicial estimado a partir del acelerómetro:", Q[0])
# Bucle para procesar todas las muestras
for t in range(1, num_samples):
    Q[t] = ekf.update(Q[t-1], gyr_data[t], acc_data[t])

print("Procesamiento completado.")
print("Forma del array de cuaterniones resultante:", Q.shape)
print("Ejemplo del primer cuaternión calculado:", Q[0])
print("Ejemplo del último cuaternión calculado:", Q[-1])