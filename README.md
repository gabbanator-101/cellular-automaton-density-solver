# cellular-automaton-density-solver
Description:
Designed and implemented a C program to simulate Elementary Cellular Automata (ECA) using custom-defined rules and solve the Density Classification Problem (DCP) through rule-based evolution. The system supports user-defined rules and initial states, dynamically simulating evolution across time steps and outputting the final classification of cell densities.
Key Features:
•	Simulated ECA using binary rules ranging from 0 to 255, applied to 3-cell neighborhoods.
•	Implemented Rule 184 (traffic flow model) and Rule 232 (majority classification) in sequence to solve the Density Classification Problem.
•	Handled array edge wrapping and binary rule interpretation for dynamic evolution of CA.
•	Developed a custom data structure to store neighborhoods and their corresponding update states.
•	Included logic to determine whether the density of ‘ON’ states is less than, greater than, or equal to 0.5 at a given time step.
•	Dynamically allocated and deallocated memory for efficient storage of evolving CA states.
•	Provided clear formatted outputs for each simulation stage to aid debugging and result interpretation.
Skills Demonstrated:
•	Bit manipulation and binary-decimal conversion.
•	Efficient memory management in C (dynamic arrays, pointers).
•	Clean code organization with modular functions and consistent naming conventions.
•	Use of assertions and input validation to ensure correctness and robustness.
