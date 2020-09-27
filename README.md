# classical-MD
Classical Molecular dynamics is simpler to code efficiently, especially on C/C++. This project involves using Green-Kubo method to compute thermal conductivity of LJ Argon.

### Details
Implementation of the Green-Kubo method is useful in computing a complex thermal property - thermal conductivity. For simplicity, the LJ potential is used on Ar liquid. The system is first equilibrated in an NVT ensemble controlled by a Nose-Hoover thermostat. The code makes the system automatically shift to an NVE ensemble for energy equilibration. This shifting is done by changing the value of eta to 0 at a particular user-defined timestep. After NVE equilibration, heat flux data is collected for every 5 timesteps. This heatflux data is post-processed to calculate the thermal conductivity and the HCACF (heat current autocorrelation function).

### Files
Main code (NVT + NVE + heatflux): NVE_NVT_thermalCond.cpp

Post processing (Thermal Conductivity): thermal_conduct.cpp

#### NOTE
Please refer to the project report (achar_Islamov.pdf) for more details regarding this project. 
