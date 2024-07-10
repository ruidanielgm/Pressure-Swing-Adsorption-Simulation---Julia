# Pressure-Swing-Adsorption-Simulation---Julia

This repository contains the code I'm using to develop a Pressure Swing Adsorption (PSA) simulation model, based on the one available here: https://www.sciencedirect.com/science/article/pii/S1383586613005091

This process consists of multiple steps. Among them are:

  1)Feed
  2)Pressurization
  3)Counter-Current Feed
  4)Purge
  5)etc.
  
The equations governing all these steps are the same and consist of mass balances for the gas and solid (adsorbent) phases, energy balances, and momentum balance (Ergun Equation). The only thing that changes between steps is the boundary conditions.

As a continuous and cyclic process, it is important to consider the variation of variables in space and time, leading to a system of PDEs, ODEs, and NAE.

The aim is to implement the code step by step:

  1)Develop the code for the Feed step (similar to a fixed bed experiment)
    1a)Implement Mass Balance
    1b)Add Momentum Balance
  1c)Add Energy Balance
  2)Develop the code for the other PSA steps
  3)Add other system particularities, including measurement delays, pipeline effects, etc.
  4)Validate the model with lab data

Completed and To-Do Tasks:

Currently, I've implemented the mass balance and momentum balance. However, it produces unexpected results when compared to a similar simulation made on gPROMS. The output is significantly delayed, and the pressure does not fall (as it should); instead, it is increasing. This seems to be the primary issue that I am looking to solve now.

If you need further assistance with specific parts of your code or have more details to share about the issues, feel free to let me know!


As an example, I'm simulating a system of two components: CH4 and CO2. The initial and the feed composition can be changed
