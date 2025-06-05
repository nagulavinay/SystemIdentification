# SystemIdentification
Investigation of Dynamic Systems Using Different Models

# Project Overview
This project focuses on identifying linear dynamic systems to form the basis for adaptive control using measured data. The main tasks involve investigating dynamic systems, estimating models through MATLAB programming, and analyzing their performance under various conditions. The project adheres to concepts from Nonlinear System Identification including system modeling, optimization, and data analysis techniques.

# Systems to Identify
**System 1:**
G1(s)= 10 / s3 + 5.8s2 + 8s + 20 (3rd-order system, PT3)

**System 2:**
G2(s)= 0.25 / s2 + 1.5s + 0.5 (2nd-order system with dead time, PT2-Tt)

**Preparatory Work in MATLAB:**
Familiarize with MATLAB programming, focusing on commands essential for system identification. The basic codes used are help(command), tf (e.g., s=tf(’s’), z=tf(’z’,T0), G=tf([b1 b0],[a2 a1 a0]), G.num, G.den, G.InputDelay, get(G), usw.), pole(G), step(G,t), impulse(G,t), lsim(G,u,t), bode(G), randn, rng, inv, fminunc, figure, plot(y,t), filter(b,a,u).

**Analyzing the given systems:**
The systems are analysed by Pole locations, step response, impulse response, and Bode diagrams. The basic codes used to analyse are s=tf(’s’), G=0.25/(s^2+1.5*s+0.5)*e^(-s), pole(G), step(G), impulse(G), bode(G)).

**Model Estimation Programs:**
Implemented MATLAB code to estimate ARX, FIR, and OE models using least squares and recursive methods.

**ARX Model (Auto Regressive with eXogeneous Input):**
Developed functions for regression matrix creation, parameter estimation, and recursive output calculation.

**FIR Model (Finite Impulse Response Model):**
Developed functions for regression matrix creation, parameter estimation, and model output calculation.

**OE Model (Output Error Model):**
Based on recursive calculation and nonlinear optimization (using fminunc) for parameter estimation.

# Model Evaluation
To assess the accuracy and robustness of the models identified for the dynamic systems, I have conducted several analyses. These evaluations focused on comparing the output of each model with the actual system output and investigating how the models perform under different conditions such as measurement noise and dead time mis-specifications.

Root Mean Square Error Calculation
Frequency Response Comparison
Effect of Measurement noise
Impact of Daed time

**Estimation of Both Systems (ARX, OE, FIR) without Measurement Noise**
Set up MATLAB transfer functions for System 1 and System 2. Generated step and impulse response data using MATLAB commands (step and impulse).

Testing without Noise: Applied the impulse and step responses to test each model’s performance, ensuring the implementation of ARX, OE, and FIR models functions correctly. Plotted the estimated model outputs alongside the actual noiseless system responses to visualize accuracy.

**APRBS-Based Estimation:**
Used Generated Input data using an APRBS (Amplitude Pseudo-Random Binary Signal) and compared each estimated model's step and impulse responses with those of the actual systems.

**Selection of Parameters:**
Appropriate sampling time T0, signal length, and model order based on system analysis. System 1 (PT3): Higher-order model, a higher model order. System 2 (PT2-Tt): Taking into Account dead time and a model order compatible with the 2nd-order behavior.

# Effect of Incorrect Dead Time in ARX Estimation for System 2
**Introduce Dead Time Error:**
Analysed ARX estimation without accurately considering the dead time or assuming an incorrect dead time (e.g., not a multiple of the sampling time T0) comparing it to the actual output to see how inaccurate dead time influences accuracy. Observed the model response and calculated RMSE to quantify the error introduced by incorrect dead time assumptions.

# Effect of Measurement Noise on ARX Estimation for System 1:
**Noise Levels and Signal Types:**
Impulse and step responses for System 1 with varying levels of measurement noise added to the output signal. Used randn function with different standard deviations (e.g., Sigma = 0.01, 0.05, 0.1) to simulate noise.

**ARX Estimation with Noise:**
Applied ARX estimation on noisy data for both step and impulse responses. Investigated how noise level, signal type, and data length affect the estimation accuracy.

# Improvement Techniques and FIR Model Behavior with Noise
**Improvement Techniques on Noisy Step Response for System 1:**

**1/A Filtering:**
Filtering the output data to reduce noise effects.

**OE Model with Nonlinear Optimization:**
Comparing with OE model for enhanced noise robustness.

**Analysis and Drawbacks:**
Examined the benefits and limitations of each technique, such as computational overhead or reduced adaptability, which technique provides the best balance between accuracy and robustness under noise.

**FIR Model with Noise:**
Test the FIR model on System 1 data with different noise levels.
Evaluated FIR’s sensitivity to measurement noise and the model's effectiveness in producing stable estimates.

**Reusable Functions:**
Developed flexible functions for creating regression matrices, estimating parameters, and calculating outputs, Ensuring the functions can accommodate varying model orders.

**RMSE Calculation:**
For each estimated model.

**Visualization**
Comparison Plots:
For each method, both the actual system outputs and the model-estimated outputs are compared.
