# Description
The code runs the results obtained for our paper [Large Deviations in Safety-Critical Systems with Probabilistic Initial Conditions](https://arxiv.org/abs/2405.13506)

This constitutes an extension for safety analysis of dynamical systems with uncertain initial conditions, so as to include unsafe events characterized by extremely low probabilities, i.e. rare events. In the presence of rare events, traditional safety methods struggle to compute the accurate probability of the event due to numerical instabilities, scalability in the system's dimension, conservativeness, etc, thereby requiring rare-event tailored methods. 

In the paper, we leverage the results from Large Deviations theory to reach an accurate and computationally feasible formulation for the probability density of initial states, which describe how probable is departing from a point and hitting an unsafe compact set. Maximizing this density allows us to determine interesting quantities, such as the most probable initial condition, hitting time and trajectory hitting the unsafe set.

The approach is demonstrated by this code using a practical example: a short-term conjunction of two space objects in quasi-circular orbit.

# Dependencies
The method relies on solving a variational problem, which we formulate in MATLAB and solve using CasADi.

You can install CasADi from [here](https://web.casadi.org/get/)
