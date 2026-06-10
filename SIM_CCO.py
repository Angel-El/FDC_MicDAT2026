# %%

import random
from scipy.signal import welch, butter, lfilter
from cco import CCO_preamp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from matplotlib.ticker import MaxNLocator
import plotly.express as px
from scipy.stats import norm

seed = 42
np.random.seed(seed)

# --- System Parameters ---
T_sim = 900e-6  # Total simulation time
Fs = 100e6
clk = 1 / Fs  # 20 ns clock period
Tint = 100e-6  # Integration time (100 µs)
T_no_int = 3.7e-6  # Non-integration time (120 µs)
clk_cco = 1.25e6
Vref = 1.25
C_cp = 270e-15
C_int = 1.72e-12
Cos = C_cp
k = 1.38e-23
T = 323.15
Rc_ota = 10e6
Cl_ota = 20e-12

# --- Variables Initialization ---
time = np.arange(0, T_sim, clk)
time_x = len(time)
# Clock cycles during active sync (Tint / clk) - 1, since the first cycle is the activation cycle.
# (Tint / clk) = 10000 clock cycles. Active integration window.
window_active_sync = round(Tint / clk)
window_gap_sync = 1

# --- Input Current ---
muestras = 300
stop = 8.425e-9
Iin = np.linspace(0, stop, muestras)
Ioffset = 8.4e-9
current_in_y = len(Iin)

# --- Noise Enables ---
enable_noise_Ioffset = 0
enable_noise_comparator = 0
enable_noise_cp = 0
enable_noise_ref = 0
enable_noise_preamp = 0
enable_noise_ota = 0

# --- Preamplifier Physical Parameters ---
I_d = 10e-6
gm_over_id = 10
gm = gm_over_id * I_d
R_d = 200e3
C_l = 200e-15
tau = R_d * C_l

# --- OTA Parameters ---
Cl_ota = 10e-12
gm_ota = 1e-3  # OTA Transconductance (1 mS)


# Function to generate the real noise curve since the psd noise of spectre analysis of each noise source
def generar_ruido_real(n_filas, n_columnas, fs, psd_white, psd_f1Hz):
    nyquist = 0.5 * fs
    freqs = np.fft.rfftfreq(n_columnas, d=1/fs)
    print(f"first medible frequency: {freqs[1]:.2f} Hz")

    # --- STEP 1: THERMAL NOISE (WHITE) ---
    thermal_raw = np.random.normal(0, 1, (n_filas, n_columnas))
    thermal_scaled = thermal_raw * np.sqrt(psd_white * nyquist)

    # =========================================================================
    # --- STEP 2: DYNAMIC AND GENERALIZED FLICKER NOISE (1/f) GENERATION -----
    # =========================================================================

    # 1. Define the geometric profile of the 1/f filter in terms of power.
    # As power decays as (1/f), the spectrum amplitude decays as 1/sqrt(f).
    with np.errstate(divide='ignore', invalid='ignore'):
        filtro_flicker = 1.0 / np.sqrt(freqs)
    # Remove the DC component (0 Hz) to avoid division by zero (infinity).
    filtro_flicker[0] = 0

    # 2. RADICAL GENERALIZATION: Find the first physical real measurable frequency.
    # It does not matter if it is 100 Hz, 1 Hz, or 0.1 Hz. It is the first useful point after DC.
    f_primera_medible = freqs[1]

    # Apply the physical law of real flicker noise (PSD = psd_f1Hz / f).
    # Calculate how much power (A²/Hz) the circuit MUST exactly have at this first frequency.
    psd_objetivo_primera_freq = psd_f1Hz / f_primera_medible

    # Evaluate how much power our baseline mathematical profile generates at that same point (bin [1]).
    # By squaring the filter amplitude, we obtain its intrinsic unscaled power.
    potencia_base_filtro_en_primer_bin = filtro_flicker[1]**2

    # The scale factor is the ratio between the target power and the default filter power.
    # The square root is used because this factor will multiply the noise AMPLITUDE.
    factor_escala = np.sqrt(psd_objetivo_primera_freq /
                            potencia_base_filtro_en_primer_bin)

    # 3. Generate the Flicker noise matrix for all rows
    flicker = np.zeros((n_filas, n_columnas))

    for i in range(n_filas):
        # Generate complex random white noise (random phase and amplitude for each frequency bin)
        blanco_en_frecuencia = np.random.normal(
            0, 1, len(freqs)) + 1j * np.random.normal(0, 1, len(freqs))

        # Shape the spectrum by multiplying it by the 1/sqrt(f) slope
        flicker_en_frecuencia = blanco_en_frecuencia * filtro_flicker

        # Inverse Fast Fourier Transform: convert the shaped spectrum back to the time domain
        flicker_temporal_puro = np.fft.irfft(
            flicker_en_frecuencia, n=n_columnas)

        # Scale the time vector to match the real physical magnitudes of current (A)
        flicker[i, :] = flicker_temporal_puro * factor_escala

    # =========================================================================
    # --- STEP 3: COHERENT SUMMATION OF THERMAL AND FLICKER NOISE ------------
    # =========================================================================

    # Combine both two-dimensional matrices (Thermal Noise + Flicker Noise)
    ruido_total = thermal_scaled + flicker

    return ruido_total


def butter_lowpass(gbw, Fs, ruido):
    nyquist = 0.5 * Fs
    b, a = butter(1, gbw / nyquist, btype='low')
    ruido_filtrado = lfilter(b, a, ruido, axis=1)
    return ruido_filtrado
# Noise of Ioffset (Flicker + White)


gbw_ota = gm_ota / (2 * np.pi * Cl_ota)
psd_white_spectre_Ioffset = 1e-27  # A^2/Hz, Value from spectre simulation
# A^2/Hz a 1 Hz, Value from spectre simulation
psd_flicker_1Hz_spectre_Ioffset = 2.9e-27
# 2.9e-27
ruido_Ioffset = generar_ruido_real(
    n_filas=Iin.shape[0],
    n_columnas=time_x,
    fs=Fs,
    psd_white=psd_white_spectre_Ioffset,
    psd_f1Hz=psd_flicker_1Hz_spectre_Ioffset,
)

Ruido_Ioffset_filtrado = butter_lowpass(
    gbw=gbw_ota,
    Fs=Fs,
    ruido=ruido_Ioffset
)

# 2. Create the input current matrix (Ipd_matrix) by adding the base input current (Iin) to the Ioffset noise.
Ipd_base = np.tile(Iin[:, np.newaxis], (1, time_x))
# Sumamos a la matriz de corriente
Ipd_matrix = Ipd_base + Ioffset + \
    (Ruido_Ioffset_filtrado * enable_noise_Ioffset)
#########

# --- Noise Sources Generation ---
# Add noise to the voltage reference (Thermal)
Vrefrms = 1e-6 * Vref
Vnref = Vrefrms * np.random.normal(size=time_x)
Vrefnoisy = Vref + (Vnref * enable_noise_ref)

# Noisy Charge Pump
Qsigma = np.sqrt(k * T * C_cp)
Vsigma = Qsigma / C_int
# Condition to enable or disable the charge pump noise
if enable_noise_cp == 1:
    Vncp = np.random.normal(0, Vsigma, size=time_x)
else:
    # Generate an array of zeros without altering randomness
    Vncp = np.zeros(time_x)

VCPnoisy = Vrefnoisy * (C_cp / C_int) + Vncp

# Create comparator noise
Vncomparator = 0.256e-3 * np.random.normal(size=time_x)
Vref_comparator_noisy = (Vrefnoisy + Vncomparator) * enable_noise_comparator


# PREAMPLIFIER NOISE
# Thermal_noise_preamp = np.sqrt(
# (16*k*T*1.15/(3*gm*4*tau))) * enable_noise_preamp


psd_white_spectre_preamp = 47e-18  # V^2/Hz, Value from spectre simulation
# V^2/Hz a 1 Hz, Value from spectre simulation
psd_flicker_1Hz_spectre_preamp = 47e-12
Noise_preamp_spectre = generar_ruido_real(
    n_filas=Iin.shape[0],
    n_columnas=time_x,
    fs=Fs,
    psd_white=psd_white_spectre_preamp,
    psd_f1Hz=psd_flicker_1Hz_spectre_preamp,
) * enable_noise_preamp

# OTA NOISE
# ota_noise = 67e-6 * enable_noise_ota
psd_white_spectre_ota = 100e-18  # V^2/Hz,Value from spectre simulation
# V^2/Hz, Value from spectre simulation
psd_flicker_1Hz_spectre_ota = 1.6e-12
ota_noise = generar_ruido_real(
    n_filas=Iin.shape[0],
    n_columnas=time_x,
    fs=Fs,
    psd_white=psd_white_spectre_ota,
    psd_f1Hz=psd_flicker_1Hz_spectre_ota,
) * enable_noise_ota

# --- Sync Generation and CCO Instantiation ---
# Create synchronization pattern (SYNC_INT_N)
pattern_sync = np.concatenate(
    [np.ones(window_active_sync), np.zeros(window_gap_sync)])
# Repeat pattern according to total simulation time
num_reps_sync = int(np.ceil(time_x / len(pattern_sync)))
patron_total_sync = np.tile(pattern_sync, num_reps_sync)[
    :time_x]  # Match time_x size
# Repeat pattern for each row of the Iin matrix
SYNC_INT_N1 = np.tile(patron_total_sync, (current_in_y, 1))

# Instantiate CCO and inject noise to simulate the CCO block
Icco_preamp = CCO_preamp()
Vcoarse_preamp_matrix, counter_preamp_matrix = Icco_preamp.integrar_fila(
    Ipd_matrix, VCPnoisy, Vref_comparator_noisy, ota_noise, Noise_preamp_spectre)


# %% 2. nCP AND nFL CALCULATION FOR EACH CORNER
# Calculate the number of times the comparator crosses the threshold within
# the integration window and the number of cycles between the first and last change.
period = window_active_sync + window_gap_sync
filas, columnas = counter_preamp_matrix.shape

# Number of complete periods (including the last incomplete one)
n_periods = int(np.ceil(columnas / period))

# Padding to make n_periods * period a multiple of period
pad_len = n_periods * period - columnas
if pad_len > 0:
    pad = np.zeros((filas, pad_len), dtype=counter_preamp_matrix.dtype)
    mat = np.hstack((counter_preamp_matrix, pad))
else:
    mat = counter_preamp_matrix.copy()

# Reshape by periods: (filas, n_periods, period)
mat_periods = mat.reshape(filas, n_periods, period)

# Only interested in the active window within each period (first window_active_sync samples)
# Shape (filas, n_periods, window_active_sync)
active = mat_periods[:, :, :window_active_sync]

# counter_CP: number of times the comparator crosses the reference threshold (1.25V) per period
counter_CP = active.sum(axis=2).astype(int)  # Shape (filas, n_periods)

# Calculate cycles within the active window (first to last inclusive)
mask = active.astype(bool)
has_any = mask.any(axis=2)  # Check if there is any 1 in the active window

# Relative indices within the active window
first_idx = np.argmax(mask, axis=2)
# Index of the last 1 (counting backward from the end of the active window)
last_idx = window_active_sync - 1 - np.argmax(mask[:, :, ::-1], axis=2)

# cycles: last - first + 1 if there are ones, else 0
counter_cycles = np.where(has_any, last_idx - first_idx + 1, 0).astype(int)

# Calculate times between the first and last CP in each integration window
# for each input current (row of the Iin_matrix)
Matriz_Tst_total = [fila * clk for fila in counter_cycles]

# %% Calculation of Converter Parameters and Results Plots
nFL = np.array(counter_cycles)
nFL_mod = nFL[:, 1:-1]  # Exclude first and last column
max_nFL = np.max(nFL_mod)
min_nFL = np.min(nFL_mod[np.nonzero(nFL_mod)])
diff_nFL = max_nFL - min_nFL
print(f"nFL Range: {diff_nFL} (min: {min_nFL}, max: {max_nFL})")

# Total activation cycles of the CP within each integration window
nCP = counter_CP
nCP_mod = nCP[:, 1:-1]
T_prima_cco = nFL_mod / (nCP_mod - 1)
Res_T_prima_cco = 1 / (nCP_mod - 1)
F_prima_cco = 1 / T_prima_cco
F_prima_cco = np.nan_to_num(F_prima_cco, nan=0.0, posinf=0.0, neginf=0.0)
columnas_validas = F_prima_cco.shape[1]
Iin_matrix_reshaped = Ipd_matrix[:, 0:columnas_validas]


# #### RESULTS AND PLOTS ###############################################################################

# %% 9. PLOT: VOLTAGE WAVEFORM AT THE OUTPUT OF THE OTA (FIRST 3 CURRENTS)
# 1. configuration of the figure with two subplots (upper for voltage waveforms, lower for SYNC signal)
# To adjust the ratio (70% analog, 30% digital)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                               gridspec_kw={'height_ratios': [2.5, 1]})

# Convertimos el tiempo a microsegundos para escala física
time_us = time * 1e6

# =========================================================================
# Upper Plot (ax1): OTA Output Voltage Waveforms
# =========================================================================
for i in range(3):
    total_current_nA = (Iin[i] + Ioffset) * 1e9
    ax1.plot(time_us, Vcoarse_preamp_matrix[i],
             linewidth=1.8,
             label=f'$I_{{in, total}}$ = {total_current_nA:.2f} nA')

ax1.set_title('Transient Voltage Waveforms and Integration Window Synchronization',
              fontsize=12, fontweight='bold')
ax1.set_ylabel('OTA Output Voltage (V)', fontsize=11)
ax1.set_ylim(1.2, Vref * 1.2)

# Grid avanzada para microelectrónica
ax1.grid(True, which='both', linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', frameon=True,
           facecolor='white', edgecolor='none')


# =========================================================================
# Downper Plot (ax2): SYNC Signal
# =========================================================================
# Extraemos el patrón SYNC correspondiente a la primera fila (o la i-ésima)
# Usamos un color profesional (p. ej., gris oscuro o azul marino) para señales digitales
ax2.step(time_us, SYNC_INT_N1[0], where='post',
         color='#2c3e50', linewidth=1.5, label='SYNC')

ax2.set_xlabel('Time ($\mu$s)', fontsize=11)
ax2.set_ylabel('SYNC (arb.u.)', fontsize=11)

# Ajustamos límites en Y para la señal digital para que no toque los bordes físicos
ax2.set_ylim(-0.2, 1.2)
ax2.grid(True, linestyle='--', alpha=0.4)


# =========================================================================
# CONFIGURACIÓN GLOBAL DE EJES Y LÍMITES
# =========================================================================
# Definimos el zoom en el eje X compartido (por ejemplo, los primeros 250 us)
plt.xlim(0, 250)

# Localizador de divisiones en el eje X
ax2.xaxis.set_major_locator(MaxNLocator(nbins=10))

plt.tight_layout()
# Ajuste fino para reducir el espacio en blanco entre las dos gráficas
plt.subplots_adjust(hspace=0.08)

plt.show()

# %% TRANSFER FUNCTION PLOT: Iin vs Fprima CCO
plt.figure(figsize=(8, 5))
plt.plot(Iin[:], F_prima_cco[:, 2], marker='.')
plt.xlabel("Input Current (A)")
plt.ylabel("D_out = F'cco = (nCP-1/nFL)")
plt.title("Transfer Function Iin VS F'cco")
plt.grid(True)
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()
np.save('xTF.npy', Ipd_base[:])
np.save('yTF', F_prima_cco[:, 2])


# %% INL CALCULATION AND PLOT FOR TRANSFER FUNCTION (Iin vs Fprima CCO)
Frec_mean = np.mean(F_prima_cco, axis=1)
Iin_mean = np.mean(Iin_matrix_reshaped, axis=1)
coefs_frec = np.polyfit(Iin_mean, Frec_mean, 1)
K_cco = coefs_frec[0]  # en Hz/A
Frec_i_fit = np.polyval(coefs_frec, Iin_mean)
INL_ppm_FSR = ((Frec_mean - Frec_i_fit) /
               (np.max(Frec_i_fit) - np.min(Frec_i_fit))) * 1e6
INL_ppm_reading = ((Frec_mean - Frec_i_fit) / Frec_mean) * 1e6
print(coefs_frec)
print(
    f"INL ppm FSR Max: {np.max(np.abs(INL_ppm_FSR)):.2f} ppm")
print(
    f"INL ppm FSR of first current: {np.abs(INL_ppm_FSR[0]):.2f} ppm")
print(
    f"INL ppm of reading Max: {np.max(np.abs(INL_ppm_reading)):.2f} ppm")
print(
    f"INL ppm of reading of first current: {np.abs(INL_ppm_reading[0]):.2f} ppm")

plt.figure(figsize=(8, 4))
plt.plot(Iin_mean, INL_ppm_FSR, marker='o',
         label='INL, ppm FSR', color='orange')
plt.xlabel("Mean Input Current (A)")
plt.ylabel("INL, ppm Full Scale Range (FSR)")
plt.title(
    "INL Curvature ")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# %% CONVERTER RESOLUTION
Delta_in = (Ipd_matrix[:, 3] / (nFL_mod[:, 3] + 1))
res_convertidor = np.log2(Ipd_matrix[:, 6] / Delta_in)
plt.plot(Iin, res_convertidor)
plt.xlabel("Input Current (A)")
plt.ylabel("Detector Resolution (bits)")
plt.title("Detector Resolution vs Input Current")
print(
    f"Converter resolution at 0nA input current: {res_convertidor[0]:.2f} bits")


# %% NOISE ANALYSIS: Qin DISTRIBUTION
# Extract the target column for noise analysis
muestras = T_prima_cco[0]
# Scale factor from Coulombs to femtoCoulombs
to_fC = 1e15
Qres = Ioffset * clk
muestras_fC = muestras * Qres * to_fC
mu = np.mean(muestras)
Qin = np.std(muestras_fC)
sigma = np.std(muestras)

# Normalized Histogram
count, bins, ignored = plt.hist(
    muestras, bins=20, density=True, alpha=0.6, color='g',
    label='Measured Data, Gaussian Fit; ' f'\n$\mu$={mu:.7f}, $Qin(fC)$={Qin:.20f}, $\sigma$={sigma:.10f}'
)

plt.xlabel('Output Code (D_out)')
plt.ylabel('Repeatability')
plt.title(f'STD Code distribution (with Ioffset, Ipd={0}), PYTHON')
plt.legend()
plt.tight_layout()
plt.show()
np.save('datos_simcco.npy', muestras)


# %% SUMMARY TABLE FOR NOISE SOURCES
# Define noise status in a structured dictionary
datos_tabla = {
    "Noise Source": [
        "Ioffset noise ",
        "OTA noise ",
        "Preamp noise ",
        "Comparator noise ",
        "Charge Pump noise ",
        "Reference noise "
    ],
    "Status (Enabled/Disabled)": [
        "Enabled" if enable_noise_Ioffset == 1 else "Disabled",
        "Enabled" if enable_noise_ota == 1 else "Disabled",
        "Enabled" if enable_noise_preamp == 1 else "Disabled",
        "Enabled" if enable_noise_comparator == 1 else "Disabled",
        "Enabled" if enable_noise_cp == 1 else "Disabled",
        "Enabled" if enable_noise_ref == 1 else "Disabled"
    ],
    "Total Noise Value (Std)": [
        f"Qin_Ioffset(fC)={Qin:.5f}" for _ in range(6)
    ]
}

# Create and display Summary DataFrame
df_resumen = pd.DataFrame(datos_tabla)
print("\n" + "="*100)
print("                    NOISE SOURCES SUMMARY (Quant. noise always present, 0.02fC)             ")
print("="*100)
print(df_resumen.to_markdown(index=False))
print("-"*100)
print(f"  Total Noise Qin(fC): {Qin:.5f}")
print("="*100)


# %% DATA EXTRAPOLATION AND FINAL VERIFICATIONS
Delta_in_exp = np.diff(Ipd_matrix, axis=0)
Tcco = nFL_mod[:, 0] * clk
Delta_nFL = np.abs(nFL_mod[1, :]-nFL_mod[0, :])
Qlsb = (Delta_in_exp[0] * Tcco[0]) / Delta_nFL[0]

print(
    f"Simulated Delta_in for the first input current: {Delta_in_exp[0, 0]:.2e} A")
print(
    f"Converter resolution at the first input current: {res_convertidor[0]:.2f} bits")
print(f"Qlsb for the first input current: {Qlsb[0] * 1e15:.4f} fC/LSB")

# %%
