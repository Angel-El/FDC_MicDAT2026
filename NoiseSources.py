# %%

from scipy.integrate import quad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
seed = 42
np.random.seed(seed)


# ==========================================
# 1. Constant Definitions
# ==========================================
time = 100000
Tobs = 10e-3  # Total observation time (10 ms)
Fs = 100e6
Tint = 100e-6  # Integration time (100 µs)
Vref = 1.25
C_cp = 270e-15
C_int = 1.72e-12
k = 1.38e-23
T = 323.15
Ioffset = 8.4e-9
Tclk = 1 / Fs

# --- OTA PHYSICAL PARAMETERS ---

# gm_ota =   # OTA Transconductance (1 mS)
R_ota = 10e6
Cl_ota = 10e-12

# --- PREAMPLIFIER PHYSICAL PARAMETERS ---

I_d = 10e-6
gm_over_id_preamp = 10
gm_preamp = gm_over_id_preamp * I_d      # Preamplifier transconductance
R_dpreamp = 200e3                # 100 kOhm
C_lpreamp = 200e-15               # 100 fF
tau_preamp = R_dpreamp * C_lpreamp
A0_preamp = gm_preamp * R_dpreamp
gamma = 2/3  # Noise factor for a MOSFET in saturation


# Noise enablers
enable_noise_quant = 1
enable_noise_Ioffset = 1
enable_noise_comparator = 1
enable_noise_cp = 1
enable_noise_ref = 1
enable_noise_preamp = 1
enable_noise_ota = 1

# Quantization Noise
Qres = Ioffset * Tclk
Q_n2_quant = (Qres/np.sqrt(12))**2
print('variance of Quantization noise, Q_n2_quant: ', Q_n2_quant)

# Charge Pump Noise
Qcp = np.sqrt(k * T * C_cp)
Vcp_rms = Qcp / C_int
Q_n2_cp = (Vcp_rms * C_int)**2
print('Variance of CP noise, Q_n2_cp: ', Q_n2_cp)


# Reference Noise
Vn_ref_rms = 1e-6*Vref
Q_n2_ref = (Vn_ref_rms * C_int)**2
print("Variance of reference noise, Q_n2_ref: ", Q_n2_ref)


# Preamplifier Noise
v_n2_preamp = (16*k*T*1.15)/(3*gm_preamp*4*tau_preamp)
Q_n2_preamp = (v_n2_preamp * C_int**2)
print("Variance of preamp noise, Q_n2_preamp: ", Q_n2_preamp)


# Comparator Noise
Vncomparator_rms = 0.256e-3
Q_n2_comp = ((Vncomparator_rms/A0_preamp) * C_int)**2
print("Variance of comparator noise, Q_n2_comp: ", Q_n2_comp)


# OTA Noise
Vn_ota_rms = 67e-6
fc_ota = 1/(2*np.pi*R_ota*Cl_ota)  # OTA corner frequency (fc)
Q_n2_ota = (Vn_ota_rms * C_int)**2
print("Variance of OTA noise, Q_n2_ota: ", Q_n2_ota)

# Ioffset Noise


def calculo_ruido_ioffset(Tint, Tobs, fc):
    # ==========================================
    # 1. Parameter definition
    # ==========================================
    S_thermal_ioffset = 1e-27  # A^2/Hz
    Kf_ioffset = 2.9e-27       # A^2/Hz at 1 Hz
    f_min = 1 / Tobs

    # ==========================================
    # 2. Integrand definition
    # ==========================================
    def integrand(f):
        # Power Spectral Density (PSD)
        S_f = S_thermal_ioffset + (Kf_ioffset / f)

        # Sinc function (avoiding division by zero at f=0 if necessary)
        if f == 0:
            sinc_val = 1.0
        else:
            sinc_val = np.sin(np.pi * f * Tint) / (np.pi * f * Tint)

        # Low-pass filter transfer function (OTA)
        lpf_ota = 1 / (1 + (f / fc)**2)

        # Filter power squared |H|^2
        H_sq = (sinc_val**2) * lpf_ota

        return S_f * H_sq

    # ==========================================
    # 3. Numerical integration to infinity
    # ==========================================
    varianza, error = quad(integrand, f_min, np.inf)
    ruido_rms = np.sqrt(varianza)

    # ==========================================
    # 4. Console output results
    # ==========================================
    print("-" * 65)
    print("         INTEGRATED RMS NOISE CALCULATION          ")
    print("-" * 65)
    print(f" Minimum frequency (f_min)  : {f_min:.3e} Hz")
    print(f" OTA Bandwidth (fc)         : {fc:.3e} Hz")
    print("-" * 55)
    print(f" Calculated Variance       : {varianza:.3e} V^2")
    print(f" Total I_offset Noise (RMS) : {ruido_rms:.3e} V")
    print("-" * 65)

    return ruido_rms


if __name__ == "__main__":
    Tint = 100e-6  # 100 us
    Tobs = 1000
    fc = fc_ota

    In_ioffset_rms = calculo_ruido_ioffset(Tint, Tobs, fc)

Q_n2_ioffset = (In_ioffset_rms * Tint)**2
print("Variance of Ioffset noise, Q_n2_ioffset: ", Q_n2_ioffset)

# Total Noise
Qtotal = np.sqrt((Q_n2_quant*enable_noise_quant) + (Q_n2_cp*enable_noise_cp) + (Q_n2_ref*enable_noise_ref) + (Q_n2_preamp*enable_noise_preamp) +
                 (Q_n2_comp*enable_noise_comparator) + (Q_n2_ota*enable_noise_ota) + (Q_n2_ioffset*enable_noise_Ioffset))

# Define individual noise states in a list
valores_individuales = [
    np.sqrt(Q_n2_quant) * 1e15 if enable_noise_quant else 0,
    np.sqrt(Q_n2_ioffset) * 1e15 if enable_noise_Ioffset else 0,
    np.sqrt(Q_n2_ota)*1e15 if enable_noise_ota else 0,
    np.sqrt(Q_n2_preamp)*1e15 if enable_noise_preamp else 0,
    np.sqrt(Q_n2_comp)*1e15 if enable_noise_comparator else 0,
    np.sqrt(Q_n2_cp)*1e15 if enable_noise_cp else 0,
    np.sqrt(Q_n2_ref)*1e15 if enable_noise_ref else 0
]
datos_tabla = {
    "Noise Source": [
        "Quantization",
        "Ioffset noise ",
        "ota noise ",
        "preamp noise ",
        "Comparator noise ",
        "Charge Pump noise ",
        "Reference noise "
    ],
    "Status (Enabled/Disabled)": [
        "Enabled" if enable_noise_quant == 1 else "Disabled",
        "Enabled" if enable_noise_Ioffset == 1 else "Disabled",
        "Enabled" if enable_noise_ota == 1 else "Disabled",
        "Enabled" if enable_noise_preamp == 1 else "Disabled",
        "Enabled" if enable_noise_comparator == 1 else "Disabled",
        "Enabled" if enable_noise_cp == 1 else "Disabled",
        "Enabled" if enable_noise_ref == 1 else "Disabled"
    ],
    "Total Noise Value (Std)": [
        f"Contribution (fC RMS): {val:.5f}" for val in valores_individuales
    ]


}

# Create and display the DataFrame
df_resumen = pd.DataFrame(datos_tabla)
# Print a clean summary
print("\n" + "="*100)
print("                                   NOISE SOURCES SUMMARY                             ")
print("="*100)
print(df_resumen.to_markdown(index=False))
print("-"*100)
print(f"  Total Noise Qin(fC): {Qtotal*1e15:.5f}")
print("="*100)


print(f"Qtotal: {Qtotal*1e15: .2f} fC")

# %%
# =====================================================================
# 2. PIECHART
# =====================================================================
# Calculate the effective contribution of each block
contribuciones_completas = {
    "Quantization": Q_n2_quant * enable_noise_quant,
    "Charge Pump": Q_n2_cp * enable_noise_cp,
    "Reference": Q_n2_ref * enable_noise_ref,
    "Preamplifier": Q_n2_preamp * enable_noise_preamp,
    "Comparator": Q_n2_comp * enable_noise_comparator,
    "OTA Integrator": Q_n2_ota * enable_noise_ota,
    "Current Offset": Q_n2_ioffset * enable_noise_Ioffset,
}

# Filter out deactivated blocks (contribution == 0) to keep the plot clean
bloques_activos = {k: v for k, v in contribuciones_completas.items() if v > 0}

labels = list(bloques_activos.keys())
valores = list(bloques_activos.values())
ruido_total = Qtotal

# =====================================================================
# 3. DRAW THE PIE CHART
# =====================================================================
plt.figure(figsize=(10, 7))

# Elegant colors for engineering reports
colores = ["#4ea8de", "#56cfe1", "#72efdd",
           "#80ffdb", "#64dfdf", "#48cae4", "#00b4d8"]

# Create the pie chart
plt.pie(
    valores,
    labels=labels,
    autopct="%1.1f%%",  # Displays percentage with one decimal place
    startangle=140,  # Rotates the chart for better readability
    colors=colores[: len(labels)],
    wedgeprops={"edgecolor": "white", "linewidth": 2, "antialiased": True},
    textprops={"fontsize": 11, "weight": "bold"},
)

# Add a clean title displaying total noise in scientific notation
plt.title(
    f"System Total Noise Contribution ($Q_{{total}}$ = {ruido_total:.3e} $C$)",
    fontsize=14,
    weight="bold",
    pad=20,
)

# Ensure the plot is rendered as a perfect circle
plt.axis("equal")

# Show the plot
plt.tight_layout()
plt.show()
