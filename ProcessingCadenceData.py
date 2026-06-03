# %% 1. DATA LOADING
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

# Input file from Cadence Spectre extraction. The Iprocessing in schematic
# print in a CSV the status of the following items
# 'SYNC_INT_N', 'out_comp', 'nFL', and 'Iin_nA' across the tran simulation.

CSV_FILE = 'datos_extraidos.csv'
df_total = pd.read_csv(CSV_FILE)
df_total.columns = df_total.columns.str.strip()

# Architecture and conversion constants
C_cp = 270e-15
Vref = 1.25
Fs = 100e6
Tint = 100e-6
gain_factor = Vref * C_cp * Fs * Tint
Tclk = 1 / Fs
Ioffset = 8.4e-9
Qres = Ioffset * Tclk

print(
    f"Data successfully loaded. Detected currents: {df_total['Iin_nA'].unique()} nA")


# %% 2. AUTOMATIC PROCESSING PER CURRENT CORNER
# This dictionary stores all extracted results grouped by current step
resultados_corners = {}

for corriente, grupo in df_total.groupby('Iin_nA'):
    grupo = grupo.reset_index(drop=True)

    # Extract dynamic signal vectors for the current corner
    sync = grupo['SYNC_INT_N'].values
    out_comp = grupo['out_comp'].values
    nfl_raw = grupo['nFL'].values

    # Integration window edge detection (SYNC_INT_N transition from 0 to 1)
    diff_s = np.diff(sync, prepend=0)
    inicios = np.where(diff_s == 1)[0]
    fines = np.where(diff_s == -1)[0]

    lista_ncp = []
    lista_nfl_efectivo = []
    # Dynamic array to store individual cycle-by-cycle nFL measurements
    lista_nfl_individuales_por_ventana = []

    for i in range(min(len(inicios), len(fines))):
        idx_s, idx_e = inicios[i], fines[i]
        oc_v = out_comp[idx_s:idx_e]
        nfl_v = nfl_raw[idx_s:idx_e]

        # --- nCP: Charge Pump Rising Edge Detection ---
        edges = np.diff(oc_v, prepend=0)
        idx_pulsos = np.where(edges > 0)[0]
        ncp_ventana = len(idx_pulsos)

        # --- nFL: Clock Cycle Extraction ---
        if ncp_ventana > 1:
            # Main macroscopic accumulation (Total active window duration)
            ciclos_efectivos = nfl_v[idx_pulsos[-1]] - nfl_v[idx_pulsos[0]]

            # Microscopic array extraction: isolated duration for every single triangle base
            nfl_individuales = np.diff(nfl_v[idx_pulsos])
        else:
            ciclos_efectivos = 0
            nfl_individuales = np.array([])

        lista_ncp.append(ncp_ventana)
        lista_nfl_efectivo.append(ciclos_efectivos)
        lista_nfl_individuales_por_ventana.append(nfl_individuales)

    # Pack metrics into the main corner database
    resultados_corners[corriente] = {
        'ncp': lista_ncp,
        'nfl_efectivo': lista_nfl_efectivo,
        # Target array containing isolated cycles
        'nfl_all': lista_nfl_individuales_por_ventana
    }
    print(
        f"✅ Processed {corriente} nA: {len(lista_ncp)} integration windows found.")


# %% 3. F_PRIMA CALCULATION PER CORNER
# F_prima_cco = (nCP - 1) / nFL
corrientes = []
f_primas = []
medias_f = []
stds_f = []
All_nfl = []

print(f"\n{'Iin (nA)':>10} | {'F_prima_cco (Mean)':>20}")
print("-" * 35)
muestras_por_corriente = {}

for corriente, datos in resultados_corners.items():
    ncps = np.array(datos['ncp'])
    nfls = np.array(datos['nfl_efectivo'])

    All_nfl.extend(nfls)
    # Safely compute f_prima for each window avoiding division by zero if nfl equals 0
    f_ventana_raw = np.divide(
        ncps - 1, nfls, out=np.zeros_like(ncps, dtype=float), where=nfls != 0)
    f_ventana = f_ventana_raw[1:-1]

    # Save complete raw sample array
    muestras_por_corriente[corriente] = f_ventana

    # Calculate statistics for the current corner
    f_prima_media = np.mean(f_ventana)
    corrientes.append(corriente)
    f_primas.append(f_prima_media)
    medias_f.append(f_prima_media)
    stds_f.append(np.std(f_ventana))

    print(f"{corriente:>8} nA | {f_prima_media:>20.6f}")

All_nfl = np.array(All_nfl)


# %% 4. DETAILED BREAKDOWN TABLE: Window and Corner Metrics
print(f"\n{'Current':>10} | {'Window':>8} | {'nCP (Pulses)':>12} | {'nFL (Cycles)':>12}")
print("-" * 55)

for corriente, datos in resultados_corners.items():
    ncps = datos['ncp']
    nfls = datos['nfl_efectivo']

    for i in range(len(ncps)):
        print(f"{corriente:>8} nA | {i:>8} | {ncps[i]:>12} | {nfls[i]:>12}")


# %% 5. PLOT: CCO TRANSFER FUNCTION
plt.figure(figsize=(8, 6))
plt.plot(corrientes, f_primas, marker='o', color='tab:blue', label="F'cco")
plt.title('CCO Transfer Function from Cadence Pre-Layout Simulation')
plt.xlabel('Input Current $I_{in}$ (nA)')
plt.ylabel("F'cco = (nCP - 1) / nFL")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()

np.save('corrientes.npy', np.array(corrientes))
np.save('f_primas.npy', np.array(f_primas))


# %% 6. PLOT: MEAN nCP PULSES VS INPUT CURRENT
ncpmean = np.mean([datos['ncp']
                  for datos in resultados_corners.values()], axis=1)

plt.figure(figsize=(10, 5))
plt.plot(corrientes, ncpmean, marker='o', color='tab:green')
plt.title('Mean nCP Pulses per Integration Window vs Input Current')
plt.xlabel('Input Current $I_{in}$ (nA)')
plt.ylabel('Mean nCP per Window')
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.show()


# %% 7. PLOT: RELATIVE STANDARD DEVIATION (%)
std_relativa = (np.array(stds_f) / np.array(medias_f)) * 100

plt.figure(figsize=(8, 4))
plt.plot(corrientes, std_relativa, 's--r')
plt.title('CCO Instability (Relative Standard Deviation)')
plt.xlabel('Input Current $I_{in}$ (nA)')
plt.ylabel('Relative Std Dev (%)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()


# %% 8. AUTOMATED INTEGRAL NON-LINEARITY (INL) EVALUATION
X_Iin = np.array(corrientes)
Y_Dout = np.array(f_primas)

# Perform linear best-fit via least squares
slope, offset = np.polyfit(X_Iin, Y_Dout, 1)
Y_ideal = slope * X_Iin + offset

# Determine absolute error and normalize against Full-Scale Range (FSR)
absolute_error = Y_Dout - Y_ideal
FSR = np.max(Y_ideal) - np.min(Y_ideal)
inl_ppm = (absolute_error / FSR) * 1e6

max_inl = np.max(np.abs(inl_ppm))
idx_max = np.argmax(np.abs(inl_ppm))

print("\n" + "="*50)
print(f"📊 LINEARITY REPORT (INL) - PRE-SILICON VERIFICATION")
print("="*50)
print(f"Ideal Gain Slope: {slope:.4e} (1/nA)")
print(f"Global Loop Offset: {offset:.4e}")
print(f"Worst-case INL:   {max_inl:.2f} ppm FSR (at {X_Iin[idx_max]:.1f} nA)")
print("="*50)

plt.figure(figsize=(8, 4))
plt.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.6)
plt.plot(X_Iin, inl_ppm, marker='s', color='darkblue',
         linewidth=2, label='Cadence Spectre INL')
plt.plot(X_Iin[idx_max], inl_ppm[idx_max], 'ro',
         label=f'Worst-case: {max_inl:.2f} ppm')
plt.title('Integral Non-Linearity (INL) Profile across Input Range')
plt.xlabel('Input Current $I_{in}$ (nA)')
plt.ylabel('INL Error (ppm FSR)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='best')
plt.tight_layout()
plt.show()

np.save('inl_cadence_ppm.npy', inl_ppm)

# %%
Iin = np.array(corrientes) * 1e-9
Delta_in_exp = np.diff(Iin, axis=0)
Tcco = All_nfl * Tclk
Delta_nFL = np.abs(All_nfl[1]-All_nfl[0])
Qlsb = (Delta_in_exp[0] * Tcco[0])/Delta_nFL


print(f"Qres of the first inpot current: {Qlsb*1e15} fC/LSB")

# %%
muestras = 1/muestras_por_corriente[0]


# GainFactor, Factor que convierte 'Output Code' a 'Input Charge' en Coulombs
# Factor para pasar de Coulombs a femtoCoulombs
to_fC = 1e15
muestras_fC = (muestras) * Qres * to_fC

mu = np.mean(muestras)
Qin_fC = np.std(muestras_fC)
sigma = np.std(muestras)


# Histograma normalizado
count, bins, ignored = plt.hist(
    muestras, bins=20, density=True, alpha=0.6, color='b', label='Measured Data, Gaussian Fit; ' f'\n$\mu$={mu:.8f}, $Qin(fC)$={Qin_fC:.8f}, $\sigma$={sigma:.8f}')

# Curva gaussiana teórica
x = np.linspace(min(muestras), max(muestras), 100)
# plt.plot(x, norm.pdf(x, mu, sigma), 'r', linewidth=2,
# label=f'Ajuste gaussiano\n$\mu$={mu:.7f}, $\sigma$={sigma:.7f}')

plt.xlabel('Output Code (D_out)')
plt.ylabel('Repetibility')
plt.title(f'STD Code distribution (with Ioffset, Ipd={0}), VAMS')
plt.legend()
plt.tight_layout()
plt.show()
np.save('datos_vams.npy', muestras)
