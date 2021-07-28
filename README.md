# Generalized-Bagnold-Model-Simulation
Generalized Bagnold Model Simulation with a Pressure-based Numerical Scheme

## Numerical methods

+ Fractional step method
+ WENO scheme
+ TVD RK3
+ THINC

## References
* Brosset L, Ghidaglia J-M, Le Tarnec L, et al. Generalized bagnold model[C]//The Twenty-third International Offshore and Polar Engineering Conference. Anchorage, Alaska: International Society of Offshore and Polar Engineers, 2013.
* Yabe T, Wang P-Y. Unified numerical procedure for compressibleand incompressible fluid[J]. Journal of the Physical Society of Japan, 1991, 60(7): 2105–2108.
* Yabe T, Aoki T. A universal solver for hyperbolic equations by cubic-polynomial interpolation i. one-dimensional solver[J]. Computer Physics Communications, 1991, 66(2–3): 219–232.
* Wei Z, Jiang Q, Nie S. A pressure-based numerical scheme for compressible–incompressible two-phase flows. Int J Numer Meth Fluids. 2021;1–16. https://doi.org/10.1002/fld.5029

## Eigen package REQUIRED!

See [Eigen](https://eigen.tuxfamily.org/) for more information about Eigen.

## Make
* make main_g for Bagnold model subjected to external body force
* make main_m for Bagnold model with an initial velocity
* make clean

## Run
* ./main_g liquid_density gas_density intial_pressure gamma scale out_file_name
  * Example: ./main_g 1000.0 1.0 1.0e5 1.4 1.0 case_1_s.ascii

* ./main_m cell_number out_file_name
  * Example: ./main_m 200 case_5_200.ascii
