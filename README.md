### Description of the function

`Hnm4lcp` is a `Matlab` function that can solve a linear complementarity
problem (LCP) in its standard form, which reads

$$
0 \leq x \perp (Mx+q) \geq 0,
$$

where $x \in \mathbb{R}^n$ is the real vector of unknowns, while $M \in
\mathbb{R}^{n\times n}$ and $q \in \mathbb{R}^n$ is the data. This
system means that the sought $x$ must be nonnegative componentwise ($x
\geq 0$ or $x_i\geq0$ for all $i\in[1:n]$), $y := Mx+q$ must be
nonnegative componentwise ($y \geq 0$ or $y_i\geq0$ for all $i\in[1:n]$)
and $x$ and $y$ must be perpendicular for the Euclidean scalar product
($x^\mathsf{T}y = 0$ or $x_1y_1+\cdots+x_ny_n=0$).

It is assumed that $M$ is nondegenerate, meaning that all its principal
minors are nonzero (i.e., $\det(M_{I,I}) \ne 0$ for all $I\subset
[1:n]$). There is no verification of this property (this is too
expensive) and there is no provision in the code to deal with a
degenerate $M$. The standard LCP has a unique solution whatever $q$ is
if and only if $M$ is a P-matrix (meaning that its principal minors are
positive: $\det(M_{I,I}) > 0$ for all $I\subset [1:n]$).

More is said on the use of the code in the introduction of
<tt>src/hnm4lcp.m</tt>, which can also be obtained by typing `help
hnm4lcp` in a `Matlab` window.

The code is assessed in section 4 of the paper (called DFG paper below)

>  Jean-Pierre Dussault, Mathieu Frappier, Jean Charles Gilbert (2024).
   'Polyhedral {Newton}-min algorithms for complementarity problems',
   Mathematical Programming (in revision). See also <a
   href="http://hal.inria.fr/hal-02306526"
   target="_blank">http://hal.inria.fr/hal-02306526</a>.

The code has been developed with the `Matlab` version `9.11.0.1837725
(R2021b) Update 2`.


### Directory content

This directory is structured as follows:
- AUTHORS: authors of the code
- README: this file
- LICENSE: license of the code
- VERSIONS: description of the successive changes in the code
- src: contains the source codes in `Matlab`,
- test: contains some test-cases for hnm4lcp and scripts to run them.


### Reproducing the paper tables

An ascii version of (part of) table 4.x (for x in [1:10]) of the <a
href="http://hal.inria.fr/hal-02306526" target="_blank">DFG paper</a>
can be obtained by running the script

>  test/certify_table(x,hnm4lcp,pathlcp,lcpsolve)

The results of the codes `Hnm4lcp` (if `hnm4lcp == true`), Pathlcp (if
`pathlcp == true`) and LCPsolve (if `lcpsolve == true`) are given, hence
depending on the values of the optional logical variables 'hnm4lcp',
'pathlcp' and 'lcpsolve'. Of course, if the results of Pathlcp and
LCPsolve are required, these pieces of software must have been
installed. At the time this text was written, these pieces of software
were available at the addresses

>  https://pages.cs.wisc.edu/~ferris/path.html<br>
>  http://github.com/erleben/num4lcp

The generation of some tables takes much computing time. For the
description of the problems and for the meaning of columns of the table,
see the <a href="http://hal.inria.fr/hal-02306526" target="_blank">DFG
paper</a>.

The problems can also be run individually using

>  `test/main.m`

but this script has not been cleaned up and is more messy.


### License

QPLicense.txt describes the QPL license under which the codes are
distributed.


### Download

The code `Hnm4lcp` can be downloaded in a number of manners:
- from personal: https://who.rocq.inria.fr/Jean-Charles.Gilbert/codes/hnm4lcp/hnm4lcp.html
- from Github:   https://github.com/gilbert-jch/hnm4lcp (this site)
- from HAL:      https://hal.science/hal-04799965v1 (soon)
- from SWH:      (soon)

------------------------------------------------------------------------

November 24, 2024
