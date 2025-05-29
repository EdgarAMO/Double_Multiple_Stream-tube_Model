# Double_Multiple_Stream-tube_Model
A double multiple stream-tube model based on the Paraschivoiu method.

The code works well for turbines with the following characteristics:
- chord-radius ratio below 0.10
- solidity lower than 0.25
- Global Reynolds 200,000 < Re < 10,000,000
- NACA0012, 0015 & 0018
- Straight or parabolic shape

The initial and final tip speed ratio have to be integers, and the step is 1. If you wish to modify this,
you must have some Fortran knowledge.

The code is not idiot-proof, if you enter wrong data the code will not work!
