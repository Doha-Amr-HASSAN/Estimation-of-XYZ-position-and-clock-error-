# Estimation-of-XYZ-position-and-clock-error-
This code performs the following:
1-	XYZ positions for all valid satellite at time 440992 and broadcast satellite clock error
2- Use the linearized GPS measurement equation to estimate the vector  dx
# ---------------------------------------------------------------------------
Two functions are concerned with input data read. First is “Broadcast_information” that reads the Keplerian parameters from the navigation file “eph.dat”. Second function is “read_measurements” which loads and read the observation file “rcvr.dat”. Then, the satellite positions are calculated based on the function “Sat_Position”.
# ---------------------------------------------------------------------------
The dx estimation is performed based on different situations by changing the elevation mask angle "cutoff_thr" and the residual threshold "res_thr". The troposphere and ionosphere corrections are not included in the estimation. Additionally, the calculation is performed using equal weights. 
Two scenarios are included, the first scenario involves all observations using large residual threshold of 1500 and a cutoff angle of 0°. The corresponding values of second scenario are 10, and 10° respectively. 
# ---------------------------------------------------------------------------
Input data:
reciever observation file: rcvr.dat
broadcast information: eph.dat
# ---------------------------------------------------------------------------
Output: it can be otanied by running "main.m"
The outputs are written to three text file format as follows: 
receiver_position.txt
satellite_position.txt
satellite_residues.txt
# ---------------------------------------------------------------------------
