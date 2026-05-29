import numpy as np


class CCO_preamp:
    """
    CCO Model with Physical Preamplifier (gm, Rd, Cl)
    """

    def __init__(self):
        # --- SYSTEM PARAMETERS ---
        self.Fs = 100e6
        self.clk = 1 / self.Fs
        self.V_init = 1.25
        self.C_int = 1.72e-12

        # Voltage reference for the preamplifier (center of its operating range)
        self.Vref = 1.25

        # --- PREAMPLIFIER PHYSICAL PARAMETERS ---
        I_d = 10e-6
        gm_over_id = 10
        self.gm = gm_over_id * I_d       # Preamplifier transconductance
        self.R_d = 200e3
        self.C_l = 200e-15
        self.R_ota = 10e6
        self.Cl_ota = 10e-12

        # Gain and Time Constant definition
        self.A0 = self.gm * self.R_d

        # Differential limits to prevent numerical noise integration limits from drifting
        self.V_diff_max = 2.7            # Represents a saturated output stage
        self.V_diff_min = -2.7

    def integrar_fila(self, Iin_matrix, Vcpnoisy, vref_comparator_noisy, ota_noise, Thermal_noise_preamp):
        """
        Executes row-by-row temporal integration of the FDC loop.

        Parameters:
        -----------
        Iin_matrix : ndarray
            Input current matrix (rows x columns).
        Vcpnoisy : ndarray
            Charge Pump reset voltage package including noise.
        vref_comparator_noisy : ndarray
            Comparator voltage noise contribution (zero-mean tracking).
        ota_noise : float
            OTA voltage noise floor (RMS value).
        Thermal_noise_preamp : float
            Preamplifier thermal noise contribution (RMS value).

        Returns:
        --------
        int_out : ndarray
            Integrated voltage profile across time (rows x columns).
        array_counter : ndarray
            Matrix tracking comparator firing events (1 if triggered, 0 otherwise).
        """
        Iin = Iin_matrix
        filas, columnas = Iin.shape

        # Discrete integration step voltage calculation
        delta_V = ((Iin / self.C_int) * self.clk)

        int_out = np.empty_like(delta_V, dtype=np.float64)
        array_counter = np.zeros_like(delta_V, dtype=np.int8)

        # Integration ramp state initialization (OTA Node)
        v = np.full(filas, self.V_init, dtype=np.float64)

        # Preamplifier differential state initialization (starts at 0V differential)
        v_preamp_diff = np.zeros(filas, dtype=np.float64)

        for c in range(columnas):
            # 1. INTEGRATION NODE VOLTAGE DROP
            v -= (delta_V[:, c])

            # 2. DIFFERENTIAL PREAMPLIFIER MODELING
            # Compute input voltage with randomized time-domain noise parameters
            v_in_preamp = (v + np.random.normal(0, ota_noise, size=filas) +
                           np.random.normal(0, Thermal_noise_preamp, size=filas)) - self.Vref

            # Target amplified output generation (Rd + gm network)
            # The signal updates through open-loop gain, while the noise floor
            # is evaluated at the output as a filtered RMS fluctuation.
            v_out_preamp_target = (self.A0 * v_in_preamp)
            v_preamp_diff = v_out_preamp_target

            # Optional clipping for algorithm numerical stability
            v_preamp_diff = np.clip(
                v_preamp_diff, self.V_diff_min, self.V_diff_max)

            # 3. THRESHOLD COMPARISON
            # Firing event triggers when the differential output crosses the noisy threshold.
            # The variable vref_comparator_noisy represents zero-mean noise fluctuations.
            comp = v_preamp_diff < (vref_comparator_noisy[c])

            # 4. CHARGE PUMP RESET COMPENSATION
            if np.any(comp):
                v[comp] += Vcpnoisy[c]
                # Upon input reset injection, the differential voltage loop tracks back towards 0V

            # Output array state storage for the current clock cycle
            int_out[:, c] = v
            array_counter[:, c] = comp.astype(np.int8)

        return int_out, array_counter

    def __call__(self, Iin_matrix, vcpnoisy, vref_comparator_noisy, ota_noise, Thermal_noise_preamp):
        return self.integrar_fila(Iin_matrix, vcpnoisy, vref_comparator_noisy, ota_noise, Thermal_noise_preamp)
