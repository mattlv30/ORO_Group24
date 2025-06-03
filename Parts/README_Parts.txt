%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

README Parts 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This folder contains a function associated with each part of the final project. Here, a brief overview of each part and the associated function is given. Each function is further explained in detail in the comments in the scripts to ensure a clearer comprehension of the code. Please note that inside each function, the input and the outputs are explained, and a brief overview of the function itself is provided.
For further details, please find in the folder “Report and assignment” a more detailed explanation about the mathematical model used. If something is not clear, feel free to reach us using the contact information in the main README file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part 1

This function computes the two impulses required for the Hohmann transfer, the Time of Flight (ToF) and the phasing time and angle, with analytical formulations. Successively, the equations of motion of the Kepler Restricted 2-Body Problem are numerically integrated. Please note that inside part1 you can find two versions of the numerical integration. The first one refers only to the integration of the elliptical transfer and is meant to assess the variation of the constant of motion during the transfer (if they stay constant, up to a reasonable error, or not). The second version performs a more extensive integration; it starts from the time t=0, then when t = phasing time, the first delta V is summed with the velocity. After t = ToF, the second impulse is considered. The final result is a numerical integration performed by steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part 2

A clarification is needed about part 2. Since in this part a random error is needed, it is computed only once, then stored in a variable and loaded into the part2 function. The function used to generate the random error and the actual values of the error are saved into the “Random error_part2” folder. rH_dott_err.mat and rH_err.mat are the two values of error used. The same values are used for different computations in order to compare different versions of the program. 
The part2 function uses the position and velocity error as initial conditions to solve the Hill-Chloessy-Whiltshire (HCW) equations. This is done first with the analytical formulation, to compute the minimum Delta V required and the Time of Flight (ToF), and then a numerical integration is used. The two results are compared. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part 3

This section computes the rendezvous with the chaser perturbed by atmospheric drag. At first, it propagates the unperturbed rendezvous in order to later compare it with the perturbed one; however, it does it in the equatorial CCS, while in part2 it was done only in the rìtarget Hill CCS. Then the perturbed rendezvous is propagated and compared to the unperturbed case. Finally, the perturbed rendezvous is propagated in the target relative CCS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part 4

In this part, a numerical model to assess the propagation of the contact forces along the manipulator is provided. Starting from the dimensions of the links and the angular displacement of the joints, the geometry of the robotic arm is defined. Then, the System Jacobian Matrix is derived. Finally, the Generalised Force Matrix is computed, and the results are printed. For further details about the mathematical models, please refer to the report in the folder “Report and assignment”. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part 5

The fifth part is composed of different scripts: pre_part5, part5_vers1, part5_vers2 and part5_vers3. The first script, pre_part5, plots graphs used for preliminary consideration of the problem. Those considerations are leveraged in part5_vers1 in order to find an impulse such that it completes the rendezvous with the chaser perturbed by drag acceleration with the same ToF as the unperturbed case.  In part5_vers2 and part5_vers3 it is used the fsolve function from the Optimisation Toolbox is used to find an impulse such that it completes the rendezvous with the chaser perturbed by drag acceleration; part5_vers2 works with fixed ToF, while part5_vers3 is a generalisation that works with various ToFs.
