# Pressure-Swing-Adsorption-Simulation---Julia

This repository contains the codes I'm using to develop a Pressure Swing Adsorption (PSA) simulation.
model available here: https://www.sciencedirect.com/science/article/pii/S1383586613005091

This process consists on multiple steps. Among them, we can find:
  1)Feed
  2)Pressurization
  3)Counter-Current Feed
  4)etc
  
The equations that govern all these steps are the same and consist on mass balance to the gas and solid (adsorbent) phase, energy balances and nmomentum balance (Ergun Equation). 
Only thing that changes between steps are the boundary conditions.

As a continuous and cyclic process, it is important to consider the variation of the variables in space and time, leading to a system of PDEs, ODEs and NAE.


The aim is to implement the code step by step:
  step 1) Develop the code for the Feed step 
    step 1a) Implement Mass Balance
    step 1b) Add Momentum Balance
    step 1c) Add Energy Balance
  step 2) Develop the code for the others PSA steps
  step 3) Add other system particularities, including measurement delays, pipelines effects, etc
  step 4) Validade the model with lab data

