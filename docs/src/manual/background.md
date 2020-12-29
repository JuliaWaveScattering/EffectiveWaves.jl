# Background

Assume that ``u(x)`` represents some field governed by a Helmholtz equation, such as a pressure field for acoustics. If ``u(x)`` is the transmitted field through a material filled with randomly placed particles, or inhomogeneities, then ``u(x)`` will be depend on the exact positions and types of particles. This makes $u(x)$ difficult to calculate and to experimentally measure exactly. For these reasons we focus on calculating the ensemble average $\langle u(x)\rangle$ over all possible configurations of the inhomogeneities.

This ensemble average can be written in a series of the form

$\langle u(x) \rangle = \sum_p \phi_p(x)$

where each $\phi_p$ satisfies a Helmholtz equation $\nabla^2 \phi_p(x) + k_p \phi_p(x) = 0$. The $k_p$ are the effective wavenumbers, and the $\phi_p$ depend on the geometry of the problem. We refer to each $\phi_p$ as a wavemode. For details see [Gower & Kristensson 2020](https://arxiv.org/pdf/2010.00934.pdf).

An important consequence of the theory, is that the wavenumbers $k_p$ depend only on the statistics and microstructure of the material. The $k_p$ do not depend on the shape of the material or the specific form of $\phi_p$. We give an example of this in [Equivalent symmetries](@ref).

## Plane waves

The simplest case is when using an incident plane wave and material in the shape of a plate or halfspace. This combination is an example of plannar symmetry (see [`PlanarSymmetry`](@ref)). In this case, $\phi_p$ has to be a plane wave. We use this case to calculate the ``k_p``, as the wavenumbers are the same for every symmetry.

In more detail, the package actually solves for the average scattering coefficient $\langle f_n\rangle (\mathbf r_1)$, of order ``n`` and for a particle centred at $\mathbf r_1$, which is then used to calculate the  average field $\langle u(x) \rangle$.   

For planar symmetry, the $\langle f_n\rangle also have to be plane waves, which allows us to  represent

$\langle f_n\rangle (\mathbf r_1) = \sum_p F_{p,n} \mathrm e^{i \mathbf k_p \cdot \mathbf r_1}$

We refer to the $F_{p,n}$ as the plane wave eigenvector of the wavenumber $k_p = \|\mathbf k_p \|$.
