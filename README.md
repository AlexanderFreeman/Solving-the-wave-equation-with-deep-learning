## Welcome to Wave equation project!

In this repository, we consider the solution of the one-dimensional normalized wave equation of the form:
$$\frac{\partial^2 u}{\partial \tau^2} - \frac{\partial^2 u}{\partial \xi^2} = 0$$

with chosen initial and boundary conditions:

$$
    \left\{
    \begin{array}{rcl}
         u (\tau,0) &=&  0,\\
         u (\tau,1) &=&  0,\\
         u (0,\xi) &=& \varphi(\xi),\\
         u_\tau(0,\xi) &=& \psi(\xi).
    \end{array}
    \right.
$$

using physics-informed neural networks (PINN).

## Acknowledgements

The work is based on the articles https://arxiv.org/abs/2006.11894 and https://arxiv.org/abs/1912.04737 and repository https://github.com/juansensio/nangs

## Copyright

Copyright 2022, Alexander Ryabov. Licensed under the Apache License, Version 2.0 (the "License"); you may not use this project's files except in compliance with the License. A copy of the License is provided in the LICENSE file in this repository.
