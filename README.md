# Brushless DC Motor Stator Optimisation 

**Sub-system 2:** Colin Laganier

The scripts are MATLAB scripts and requires the Global Optimization Toolbox. The execution time is for a 2.9 GHz Dual-Core Intel Core i5 (6267U) running MacOS 11.0.1. 

These scripts optimise the torque to weight ratio of brushless DC motor's stator based on 5 variables: number of rotor poles, number of stator teeth, number of coil windings, stator radius and stator thickness.

## Files:

- stator_info.m : script to output the torque to weight ratio of the stator depending on the variable input
- stator_fmincon.m :  optimisation script using fmincon with three different algorithms (interior-point, active set & SQP)
  - execution time : 0.5-1 second
- stator_geneticAlgorithm.m : optimisation script using a genetic algorithm. The returned parameters change with every run due to the random inital values selected by the algorithm.
  - execution time : 2 seconds
- stator_gamultiobj.m : optimisation script using a multiobjective genetic algorithm (objective 1: torque, objective 2: mass). Outputs the Pareto Front plot. Requires to run the second cell to select the final optimised parameters. The returned parameters change with every run due to the random inital values selected by the algorithm.
  - execution time : 20 seconds
- stator_patternSearch.m : optimisation script using a multiobjective pattern search algorithm (objective 1: torque, objective 2: mass). Uses the fmincon optimised parameters of the individual functions as initial points for the algorithm. Outputs the Pareto Front plot. Requires to run the second cell to select the final optimised parameters. 
  - execution time : 11 seconds
