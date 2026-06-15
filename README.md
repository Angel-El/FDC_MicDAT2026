# FDC_MicDAT2026
Python behavioral engine for top-down verification of an FDC loop for CT readouts.Python behavioral engine for top-down verification of an FDC loop for CT readouts.
# FDC-CT-Readout-Behavioral-Model

This repository contains the custom, clock-cycle-accurate Python behavioral simulation framework developed for the pre-silicon verification of a charge-balancing Frequency-to-Digital Converter (FDC) loop. 

This framework accompanies the paper submitted to **MicDAT 2026**: *"Analysis and Pre-Silicon Verification of a FDC based CT Detector IC With Integrated Photodiode"*.


## Repository Script Architecture

To facilitate the pre-silicon verification and post-layout benchmarking of the FDC loop, the codebase is structured into the following specialized Python modules:

*   **`SIM_CCO.py`**: The main execution script of the behavioral simulation framework. It configures the global loop setup, controls the transient simulation execution, sweeps the input current parameters, and evaluates the top-level performance metrics of the converter.
*   **`cco.py`**: A core object-oriented library that defines the internal behavioral blocks of the Current-Controlled Oscillator (CCO) and the feedback network. It explicitly models the clock-cycle-accurate dynamics of the sub-blocks, including charge integration, active amplification, high-speed comparison, and charge-pump (CP) synchronous current injection. This file is automatically imported and executed by `SIM_CCO.py`.
*   **`NoiseSources.py`**: A theoretical analytical tool dedicated to calculating the mathematical noise floor contributions of individual analog frontend blocks. It estimates noise power spectral densities to predict architectural performance limiters.
*   **`ProcessingCadenceData.py`**: A post-processing utility designed to parse and evaluate transient simulation data exported as a `.csv` matrix from **Cadence Spectre**. By importing raw foundry data (input photodiode currents, $n_{CP}$ pulse counts, $SYNC$ signal states, and transient time steps), this script extracts the fine frequency-to-digital metrics ($n_{FL}$, time between first and last nCP within integration window). This enables a direct, high-fidelity accuracy benchmarking between the theoretical Python behavioral model and the physical transistor-level implementation.


## 🚀 Key Features
* **Clock-Cycle-Accurate Simulation:** Discretizes continuous-time integration into synchronous temporal steps ($T_{\text{clk}} = 10\text{ ns}$), completely bypassing matrix-based differential equations.
* **Encapsulated Circuit Dynamics:** Internalizes CCO node dynamics, preamplifier open-loop gain modeling, and voltage-domain noise injection within a single object-oriented class execution.
* **Automated Verification Flow:** Seamless parsing and post-processing engine optimized to handle both behavioral internal tracking and external Spectre transient simulation logfiles. The complete execution pipeline matches the flow illustrated in **Fig. 2** of the paper.

---
The entire simulation environment is driven by a matrix execution approach structured around two main configuration levers: the **Total Simulation Time (T_sim)** and the **Number of Input Currents (M)**. 

Since the external synchronization window (SYNC) is fixed at T_int = 100 us, simulating a total time of T_sim = 1 ms will automatically yield 10 consecutive integration conversion cycles. Furthermore, given the main clock frequency of 100 MHz, each individual 100 us integration window precisely evaluates 10,000 clock cycles (N samples). By tweaking these parameters, the testbench generates an input data matrix of size M x N—where M represents the input current steps and N represents the total discrete temporal samples. Adjusting this matrix geometry allows you to compute all the target metrics:


## 📊 How to Reproduce Paper Python Model Results

The simulation framework is highly parameterized. To replicate the exact quantitative verification benchmarks described in the paper, adjust the top-level testbench ("SIM_CCO.py") script by modifying the total simulation time (T_sim) and the number of input currents (and the step) while keeping the physical integration window fixed at T_int = 100 us:

### 1. Intrinsic INL (Integral Non-Linearity) Evaluation
* **Objective:** Capture structural loop non-linearity across the entire dynamic range (0 nA to 175 nA).
* **Recommended Settings:** First, you have to disable all the noise sources. Then, configure a fine input current sweep array (e.g., set 500 current points (muestras)).
* **Total Simulation Time (T_sim):** Set T_sim = 1 ms (which corresponds to simulating 10 consecutive integration windows). This provides enough statistical averaging for macroscopic structural metrics while keeping the sweep computationally efficient.

### 2. Charge Resolution (Qres / Qlsb) Verification
* **Objective:** Validate the precise quantization packet sizing and threshold-crossing dynamics.
* **Recommended Settings:** Maintain the highest time-step precision configuration to minimize discrete quantization errors during comparator firing evaluations. set stop at 8.425 nA and 300 currents ("muestras = 300"), to have a fine delta input current.
* **Total Simulation Time (T_sim):** A long simulation run is **not required**. Set T_sim = 900 us (simulating just the very first integration window). The framework extracts the exact charge resolution metrics directly from this initial window execution.

### 3. Input-Referred Noise Floor Analysis
* **Objective:** Characterize worst-case charge noise under dark conditions (Ipd = 0 nA).
* **Recommended Settings:** First, enable all input noise sources. Then, configure a coarse array samples (e.g., 10 input current steps), since the noise floor analysis focuses exclusively on the first operational point fixed at the hardware baseline pedestal current (Itotal = Ioffset = 8.4 nA), worst case.
* **Total Simulation Time (T_sim):** Maximize the simulation run to T_sim = 10 ms or more (simulating 100+ consecutive integration windows). Accumulating a longer temporal profile is crucial to achieve the statistical realism required for accurate time-domain noise abstraction.

## 🔌 How to Reproduce Paper Spectre Results

The transistor-level verification data published in the paper was obtained via transient simulations in Cadence Spectre. During the simulation run, the sessional transient raw waveforms were recorded into the Spectre simulation log, which was subsequently parsed and exported into standalone CSV datasets. 

In the repository, you will find three distinct CSV files—one for each specific circuit benchmark discussed in the paper. To process and visualize these Cadence results, you must use the **`ProcessingCadenceData.py`** script. Inside this script, you can selectively toggle which database to read by modifying the target file string variable (e.g., updating the script to point to `CSV_FILE_noise`, `CSV_FILE_Chargeresolution`, or `CSV_FILE_INL`). Running the script with the corresponding file will automatically compute the data and reproduce the exact transistor-level figures and quantitative metrics presented in the paper.

### 1. Intrinsic INL (Integral Non-Linearity) Evaluation
* **Objective:** Verify structural loop non-linearity and transfer curve accuracy at the physical transistor level across the dynamic range.
* **Spectre Simulation Setup:** A parametric transient sweep was conducted utilizing an array of 20 input currents. For this particular verification metric, the absolute number of consecutive integration windows per current is not the critical factor; instead, the priority is maintaining a high transient resolution with dense time-stepping and maximum analytical points during the conversion.
* **Target Script Variable:** Set the input path in the processing script to point to `CSV_FILE_INL`.

### 2. Charge Resolution Verification
* **Objective:** Demonstrate the minimum detectable charge step and validate the sub-nanoampere resolving capabilities of the CCO architecture.
* **Spectre Simulation Setup:** A dedicated transient simulation was executed using a short integration time span (a reduced number of total integration windows). To strictly evaluate the physical resolution floor, the input stimulus was configured with extremely micro-scaled current steps separated by a delta in the picoampere range, allowing clear observation of the quantization and switching thresholds.
* **Target Script Variable:** Set the input path in the processing script to point to `CSV_FILE_Chargeresolution`.

### 3. Noise Floor Assessment (Worst-Case at 0 nA)
* **Objective:** Quantify the total charge noise floor and temporal jitter performance under dark-current or zero-signal baseline conditions.
* **Spectre Simulation Setup:** A single continuous transient simulation was executed strictly at 0 nA input current, which represents the lowest signal-to-noise ratio and the worst-case operational scenario for noise injection. To achieve adequate statistical averaging and reliable noise spectral density parsing, the simulation length was extended to a substantial duration of 10 consecutive integration windows (Tsim = 1 ms).
* **Target Script Variable:** Set the input path in the processing script to point to `CSV_FILE_noise`.

***

> ⚠️ **P.S. (Important Note on Script Execution):** Please note that since all benchmarks share the same unified post-processing framework (`ProcessingCadenceData.py`), attempting to plot specific multi-point curves while analyzing a single-point dataset will trigger dimension or interpolation errors in Python. For instance, when evaluating the **Noise Floor (`CSV_FILE_noise`)**, only a single static current (0 nA) is available; consequently, execution of code blocks intended for sweeping functions—such as the INL curve or transfer characteristics—must be skipped, as they inherently require multiple current steps. When replicating results, please focus exclusively on the specific performance metrics and figures designated for each independent simulation setup.




## 💻 Prerequisites & Installation

To run the simulation engine, ensure you have Python 3.x installed along with the following standard scientific computing libraries:

```bash
pip install numpy matplotlib

