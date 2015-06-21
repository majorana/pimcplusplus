For a particle in a box with periodic boundary conditions, the quantized energy levels in a cube of side  are

E_{ijk} = \frac{4 \hbar<sup>2 \pi</sup>2 (i<sup>2 + j</sup>2 + k<sup>2)}{2mL} = \lambda \frac{4 \pi</sup>2 (i<sup>2 + j</sup>2 + k^2)}{L}_

where

$\lambda=\frac{\hbar^2}{2m}$

The expectation value of the energy is
$\frac{1}{Z}\sum_{ijk} \left< \Psi_{ijk}|\hat H \exp[-\beta \hat H]|\Psi_{ijk} \right> = \frac{1}{Z}\sum_{ijk} E_{ijk}\exp[-\beta E_{ijk}]$
where
$Z=\sum_{ijk} \left< \Psi_{ijk}|\exp[-\beta E_{ijk}]|\Psi_{ijk} \right>$

This sum will converge rapidly to  and can be computed numerically with a few lines of Python code, for example.