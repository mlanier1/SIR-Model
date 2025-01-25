# SIR-Model

Simulation of the SIR and SEIR Models using Euler, RK2, and RK4 methods for analyzing disease propagation during the Hong Kong Flu pandemic.

## Project Overview
This project simulates the spread of the Hong Kong Flu (H3N2) pandemic using the following:
- SIR (Susceptible-Infectious-Recovered) model
- SEIR (Susceptible-Exposed-Infectious-Recovered) model

The scripts compute disease progression numerically using:
1. **Euler's Method**  
2. **2nd Order Runge-Kutta (RK2)**  
3. **4th Order Runge-Kutta (RK4)**  

The outputs include numerical approximations and insights into disease behavior under various conditions.

## Files in Repository
- **euler_sir_model.f90**: Implements the Euler method for solving the SIR model equations.
- **rk2_sir_model.f90**: Implements the 2nd Order Runge-Kutta (RK2) method for solving the SIR model equations.
- **rk4_sir_model.f90**: Implements the 4th Order Runge-Kutta (RK4) method for solving the SIR model equations.

## How to Compile and Run
1. Use a Fortran compiler such as `gfortran`.
2. To compile any of the scripts:
   ```bash
   gfortran euler_sir_model.f90 -o euler_sir
   ./euler_sir
