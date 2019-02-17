---
layout: post
title: Quantum Monte Carlo Primer
---
{% include mathjax.html %}

This post is the first in a series of blog posts detailing my experiences developing a quantum Monte Carlo library in Rust. This first entry serves as an introduction to quantum Monte Carlo methods, and a little bit of quantum mechanics itself, hopefully in a way that is understandable to the average programmer.

## Quantum Mechanics

### The Basics

The objective of quantum mechanics is to solve the Schrödinger equation (SE):

$$
  \hat{H}\Psi(\mathbf{x}, t) = i\hbar\frac{\partial \Psi(\mathbf{x}, t)}{\partial t}.
$$

On the left-hand side, we have the Hamiltonian operator \\(\hat{H}\\), acting on the *wave function* \\(\Psi\\). On the right-hand side, we have the time evolution of the wave function. All you really need to know is that $\hat{H}$ is a linear operator (in the sense of \\(d/dx\\), or a matrix), which encodes the energetic configuration of the system under consideration. For example, if we're looking at the quantum behavior of a magnet, \\(\hat{H}\\) will somehow involve the spins of each particle. In a molecular system - which is what we'll be mostly concerned with - it will involve the positions of the electrons and nuclei that make up the molecule. 

All in all, the Schrödinger equation tells us that the time evolution of quantum mechanical particle is governed by its energy. This is actually exactly the case in classical mechanics, although the equations are a little different. One important difference is that, where the classical state of a system is encoded by the positions and momenta \\(\mathbf{x}_i\\) and \\(\mathbf{p}_i\\) of the particles in the system, the quantum state is fully encoded by the wave function \\(\Psi\\). What \\(\Psi\\) exactly means is not important; all that we need to know is that given \\(\Psi\\), we have all the information about the system that we possibly can.

Now, how do we go about solving for \\(\Psi\\)? The first step is to separate the wave function into a time-dependent part and a time-independent part. I'll skip the details, but the result is that the full solution to the SE is:

$$
  \Psi(\mathbf{x}, t) = \sum_i c_i \psi_i(\mathbf{x})e^{-iE_i t/\hbar}.
$$

Here, the functions \\(\psi_i\\) and the numbers \\(E_i\\) are solutions to the *time-independent Schrödinger equation* (TISE):

$$
  \hat{H}\psi_i(\mathbf{x}) = E_i\psi_i(\mathbf{x}).
$$

The solutions to the TISE are called the *eigenfunctions* or *eigenstates* of \\(\hat{H}\\), and \\(E_i\\) are the energy eigenvalues of these states. If you've ever taken a linear algebra course, this should be familiar territory: the TISE is nothing more than an eigenvalue problem! For some more terminology, the function \\(\psi_i\\) with the lowest associated \\(E_i\\) is called the *ground state* of the system.

The problem of solving the SE then reduces to finding the solutions to the TISE, finding the expansion coefficients \\(c_i\\) (for which there is a standard recipe), and forming their linear combination with the appropriate time evolution factors \\(\exp(-E_i t/\hbar)\\). Easy, right? Well...

### The Electronic Structure Problem

While the solution to the TISE is easy to find for some very simple systems, unfortunately real life is messy and complicated, and analytical solutions are very much impossible in the most useful cases. In principle, though, we have a recipe for solving any quantum mechanical system, so let's see how bad it gets. As I said, we'll mostly look at molecular simulations, so what does the TISE look like for a molecule containing \\(N_n\\) nuclei with charge \\(Z_i\\), and \\(N_e\\) electrons? For simplicity, I'll assume that the nuclei are stationary, which is a fair approximation considering that they are far more massive than the electrons. The TISE then comes down to:

$$
  -\frac{1}{2}\sum_{i=1}^{N_e} \nabla_i^2 \psi(\mathbf{x}) - \sum_{i=1}^{N_e}\sum_{j=1}^{N_n}\frac{Z_j}{|\mathbf{R}_j - \mathbf{x}_i|}\psi(\mathbf{x}) + \sum_{i=1}^{N_e}\sum_{j > i}^{N_e}\frac{1}{|\mathbf{x}_i - \mathbf{x}_j|}\psi(\mathbf{x}) = E\psi(\mathbf{x}).
$$

Gulp.

So I don't think I need to explain why this equation is far too hopelessly complicated to solve analytically. Even with prescribed nuclear positions \\(\mathbf{R}_i\\), there is little hope of obtaining an exact, closed-form solution. Instead, we have to turn to clever approximate solution methods.

"Aha!", I hear you think. "I've taken a course on numerical solutions to PDEs! All we have to do is set up a spatial grid in \\(3N_e\\) dimensions, figure out some implicit scheme, and we're done!". If only life were so simple. If we want to get a certain prescribed accuracy, we need to use \\(n\\) grid points for each spatial dimension. Since there are \\(3N_e\\) dimensions, it's easy to see that a simple PDE solving scheme nets us an exponential time complexity in $N_e$, any computational scientist's worst nightmare. We'll have to be a bit smarter than that.

### The Variational Principle

The most useful eigenstate of \\(\hat{H}\\) is generally its ground state \\(\psi_1\\). Most molecules in nature will be found in approximately this state, meaning that the expansion coefficient \\(c_1 > c_i\\), \\(i = 2, 3, \dots\\). That cuts down our work significantly, but we still have this exponential complexity issue. As it turns out, while we can't find the numerically exact solution by directly solving the TISE, we can instead "guess" a solution. There then exists a theorem that tells us how "good" our solution is:

> For any trial function \\(\psi_T\\) that conforms to the boundary conditions imposed by our system, we have that
>
>$$
>  \langle H \rangle_{\psi_T} \equiv \frac{\int d^{3N_e}\mathbf{x} \psi_T^*(\mathbf{x}) \hat{H}\psi(\mathbf{x})}{\int d^{3N_e}\mathbf{x}|\psi_T(\mathbf{x})|^2} \geq E_1,
>$$
>
> where the inequality is satisfied only if \\(\psi_T = \psi_1\\). 

This theorem is called the variational principle. So while we can't find \\(\psi_1\\) directly, we can make an educated guess, and compute the above quantity. We can then attempt to improve our guess, and see if we reduced the value of \\(\langle H \rangle\\) for this new \\(\psi_T\\). Then we can continue this process until we see no improvements, and declare our trial function locally optimal. 

This is a lot better: instead of solving a \\(3N_e\\)-dimensional PDE, we just have to compute two \\(3N_e\\)-dimensional integrals iteratively. Unfortunately, this is still plagued by the same problem of dimensionality: if we use a grid-based numerical integraion algorithm, like Simpson's rule, the accuracy falls off with the number of dimensions in - again - an exponential fashion. Luckily, there exist methods for computing high-dimensional integrals whose accuracy does not depend on the number of dimenions: Monte Carlo integration.
