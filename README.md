# SystemIdentification
Investigation of Dynamic Systems Using Different Models

# Project Overview

This project focuses on identifying linear dynamic systems to form the basis for adaptive control using measured data. The main tasks involve investigating two given dynamic systems, estimating models through MATLAB programming, and analyzing their performance under various conditions. The project adheres to concepts from Nonlinear System Identification including system modeling, optimization, and data analysis techniques.
# Systems to Identify
# System 1: 
  G1(s)= 10 / s3 + 5.8s2 + 8s + 20 (3rd-order system, PT3)
# System 2: 
  G2(s)= 0.25 / s2 + 1.5s + 0.5 (2nd-order system with dead time, PT2-Tt)
# Preparatory Work in MATLAB
Familiarize with MATLAB programming, focusing on commands essential for system identification.

# Briefly analyze the given systems:
Pole locations, step response, impulse response, Bode diagrams.
Create MATLAB functions to generate datasets (input-output signals) for system identification.

# Model Estimation Programs
Implement MATLAB programs to estimate ARX, FIR, and OE models using least squares and recursive methods.
# ARX Model (Auto Regressive with eXogeneous Input): 
Develop functions for regression matrix creation, parameter estimation, and recursive output calculation.
# FIR Model (Finite Impulse Response Model): 
Implement functions for regression matrix creation, parameter estimation, and model output calculation.
# OE Model (Output Error Model): 
Use recursive calculation and nonlinear optimization (using fminunc) for parameter estimation.

# Model Evaluation
Evaluate model quality by comparing the model output with actual system output.
Calculate RMSE (Root Mean Square Error) for assessment.
Compare frequency responses of the real system and the model, and examine effects of measurement noise and dead time.

# 1. Estimation of Both Systems (ARX, OE, FIR) without Measurement Noise
Set up MATLAB transfer functions for System 1 and System 2.
Generate step and impulse response data using MATLAB commands (step and impulse).

Testing without Noise:
Apply the impulse and step responses to test each model’s performance, ensuring the implementation of ARX, OE, and FIR models is functioning correctly.
Plot the estimated model outputs alongside the actual noiseless system responses to visualize accuracy.
# APRBS-Based Estimation:

Generate input data using an APRBS (Amplitude Pseudo-Random Binary Signal).
Use the APRBS input data for model identification.
Compare each estimated model's step and impulse responses with those of the actual systems.
# Selection of Parameters:

Choose appropriate sampling time T0, signal length, and model order based on system analysis:
System 1 (PT3): Higher-order model, so select a higher model order.
System 2 (PT2-Tt): Account for dead time and use a model order compatible with the 2nd-order behavior.

# 2. Effect of Incorrect Dead Time in ARX Estimation for System 2
# Introduce Dead Time Error:

Use ARX estimation without accurately considering the dead time or assuming an incorrect dead time (e.g., not a multiple of the sampling time T0).
Examine the ARX model output under these conditions, comparing it to the actual output to see how inaccurate dead time influences accuracy.
Observe the model response and calculate RMSE to quantify the error introduced by incorrect dead time assumptions.
# 3. Effect of Measurement Noise on ARX Estimation for System 1:
# Noise Levels and Signal Types:

Generate impulse and step responses for System 1 with varying levels of measurement noise added to the output signal.
Use randn with different standard deviations (e.g., Sigma = 0.01, 0.05, 0.1) to simulate noise.
# ARX Estimation with Noise:

Apply ARX estimation on noisy data for both step and impulse responses.
Investigate how noise level, signal type, and data length affect the estimation accuracy.
Quantify estimation errors using RMSE and plot the responses to visualize deviations from the noiseless model.
# 4. Improvement Techniques and FIR Model Behavior with Noise
Improvement Techniques on Noisy Step Response for System 1:

# 1/A Filtering: 
Filter the output data to reduce noise effects.
# Instrumental Variables (IV) Method: 
Implement IV to handle correlation between noise and input.
# OE Model with Nonlinear Optimization: 
Compare with OE model for enhanced noise robustness.
# Analysis and Drawbacks:

Examine the benefits and limitations of each technique, such as computational overhead or reduced adaptability.
Determine which technique provides the best balance between accuracy and robustness under noise.
FIR Model with Noise:

# Test the FIR model on System 1 data with different noise levels.
Evaluate FIR’s sensitivity to measurement noise and the model's effectiveness in producing stable estimates.
Determine a suitable FIR model order for both systems based on error analysis and performance.

# Reusable Functions: 
Develop flexible functions for creating regression matrices, estimating parameters, and calculating outputs. Ensure the functions can accommodate varying model orders.
# RMSE Calculation: 
For each estimated model. 

# Visualization
# Comparison Plots: 
For each method, plot both the actual system outputs and the model-estimated outputs for easy comparison.
# Frequency Response Analysis: 
Use Bode plots to compare model frequency responses with those of the real system.
