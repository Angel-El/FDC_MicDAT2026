# FDC_MicDAT2026
Python behavioral engine for top-down verification of an FDC loop for CT readouts.Python behavioral engine for top-down verification of an FDC loop for CT readouts.
# FDC-CT-Readout-Behavioral-Model

This repository contains the custom, clock-cycle-accurate Python behavioral simulation framework developed for the pre-silicon verification of a charge-balancing Frequency-to-Digital Converter (FDC) loop. 

This framework accompanies the paper submitted to **MicDAT 2026**: *"Analysis and Pre-Silicon Verification of a FDC based CT Detector IC With Integrated Photodiode"*.

## 🚀 Key Features
* **Clock-Cycle-Accurate Simulation:** Discretizes continuous-time integration into synchronous temporal steps ($T_{\text{clk}} = 10\text{ ns}$), completely bypassing matrix-based differential equations.
* **Encapsulated Circuit Dynamics:** Internalizes CCO node dynamics, preamplifier open-loop gain modeling, and voltage-domain noise injection within a single object-oriented class execution.
* **Automated Verification Flow:** Seamless parsing and post-processing engine optimized to handle both behavioral internal tracking and external Spectre transient simulation logfiles. The complete execution pipeline matches the flow illustrated in **Fig. 2** of the paper.

---
The entire simulation environment is driven by a matrix execution approach structured around two main configuration levers: the **Total Simulation Time (T_sim)** and the **Number of Input Currents (M)**. 

Since the external synchronization window (SYNC) is fixed at T_int = 100 us, simulating a total time of T_sim = 1 ms will automatically yield 10 consecutive integration conversion cycles. Furthermore, given the main clock frequency of 100 MHz, each individual 100 us integration window precisely evaluates 10,000 clock cycles (N samples). By tweaking these parameters, the testbench generates an input data matrix of size M x N—where M represents the input current steps and N represents the total discrete temporal samples. Adjusting this matrix geometry allows you to compute all the target metrics:


## 📊 How to Reproduce Paper Results

The simulation framework is highly parameterized. To replicate the exact quantitative verification benchmarks described in the paper, adjust the top-level testbench script by modifying the total simulation time (T_sim) while keeping the physical integration window fixed at T_int = 100 us:

### 1. Intrinsic INL (Integral Non-Linearity) Evaluation
* **Objective:** Capture structural loop non-linearity across the entire dynamic range (0 nA to 175 nA).
* **Recommended Settings:** Configure a fine input current sweep array (e.g., set 500 current points (muestras)).
* **Total Simulation Time (T_sim):** Set T_sim = 1 ms (which corresponds to simulating 10 consecutive integration windows). This provides enough statistical averaging for macroscopic structural metrics while keeping the sweep computationally efficient.

### 2. Charge Resolution (Qres / Qlsb) Verification
* **Objective:** Validate the precise quantization packet sizing and threshold-crossing dynamics.
* **Recommended Settings:** Maintain the highest time-step precision configuration to minimize discrete quantization errors during comparator firing evaluations. set stop at 9nA and 60 currents, to have a fine delta input current.
* **Total Simulation Time (T_sim):** A long simulation run is **not required**. Set T_sim = 700 us (simulating just the very first integration window). The framework extracts the exact charge resolution metrics directly from this initial window execution.

### 3. Input-Referred Noise Floor Analysis
* **Objective:** Characterize worst-case charge noise under dark conditions (Ipd = 0 nA).
* **Recommended Settings:** Firs, enable all input noise sources. Then, configure a coarse array of just a few current samples (e.g., 5 input current steps), since the noise floor analysis focuses exclusively on the first operational point fixed at the hardware baseline pedestal current (Itotal = Ioffset = 8.4 nA), worst case.
* **Total Simulation Time (T_sim):** Maximize the simulation run to T_sim = 6 ms or more (simulating 60+ consecutive integration windows). Accumulating a longer temporal profile is crucial to achieve the statistical realism required for accurate time-domain noise abstraction.

## 💻 Prerequisites & Installation

To run the simulation engine, ensure you have Python 3.x installed along with the following standard scientific computing libraries:

```bash
pip install numpy matplotlib

