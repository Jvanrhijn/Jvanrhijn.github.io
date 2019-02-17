---
layout: post
title: Quantum Monte Carlo Primer part 2
---
{% include mathjax.html %}

This post is the second in a series of blog posts detailing my experiences developing a quantum Monte Carlo library in Rust. This first entry serves as an introduction to quantum Monte Carlo methods, and a little bit of quantum mechanics itself, hopefully in a way that is understandable to the average programmer.

## Quantum Monte Carlo

### Monte Carlo Integration

In the previous post, we outlined a method for finding solutions to the time-independent
Schrödinger equation. Now we will address the issue of dimensionality, using Monte Carlo
integration methods.

First, recall the formula for the average value of a function on an interval 
\\(I = \[a, b\]\\):

$$
  \langle f \rangle_I = \frac{1}{b - a}\int_a^b f(x)dx.
$$

We can turn this equation, and use an estimate of the average value of \\(f\\) to compute
the integral. In it's simplest form, we uniformly sample random points from \\(I\\), 
and compute the average of these points, yielding:

$$
  \langle f \rangle_I \approx \bar{f} = \frac{1}{N}\sum_{i=1}^N f(x_i).
$$

Given enough points, we will see \\(\bar{f} \rightarrow \langle f \rangle_I\\). What's
more, we have that the error in this approximation converges as \\( 1/\sqrt{N}\\). Now,
this exact same process holds for higher-dimensional integrals, and the real kicker is
that the error converges at the same rate, **irrespective of the dimensionality of the
integral!** 

Now, there is one extra complication here. If we sample points \\(x_i\\) uniformly, the
integral may converge very slowly if \\(f\\) is zero over a large portion of \\(I\\).
To be more efficient, we should instead sample points from a probability distribution
that that converges to \\(f\\). So we assume some distribution \\(\rho(x)\)), and write

$$
  (b - a)\langle f \rangle_I = \int_I \rho(x) \frac{f(x)}{\rho(x)}dx.
$$

If we now sample \\(f(x)/\rho(x)\\) according to the distribution \\(\rho\)), the value
of the integral clearly remains the same, but it will be sampled much more efficiently.

In practice, we sample points from \\(I\\) using a Markov chain.
Given a point \\(x_k\\), we propose an
\\(x_{k+1}\\). We then either accept this point as the next point in the chain or reject
it, according to a probability \\( P(x_{k+1}|x_k) \\). There are some requirements on 
\\(P\\) that I won't go into; suffice to say that we can factorize \\(P\\) into a trial
probability \\(T\\) and an acceptance probability \\(A\\). Given a trial probability, 
which depends on the method by which points are proposed, we can then compute the
acceptance as

$$
  A(x_{k+1}|x_k) = \min\left\{ 1, \frac{T(x_{k}|x_{k+1})\rho(x_{k+1})}{T(x_{k+1}|x_k)\rho(x_{k})} \right\}.
$$

We then generate a random number \\(r \in \[0, 1)\\). If \\(A > r\\), we accept the
proposed \\(x_{k+1}\\); otherwise, we reject it, and propose a new point. This algorithm
is called the *Metropolis-Hastings algorithm*. 

While I used one-dimensional integration as an example, the entire algorithm holds just
as well for multidimensional integrals.

### Monte Carlo for Variational Integrals

Now that we have a recipe for computing multidimensional integrals accurately, let's
apply it to the problem at hand. Recall the variational integral,

$$
  \langle H \rangle = \frac{\int \psi^* (x)\hat{H}\psi(x) d^{3N}x}{\int \psi^* (x)\psi(x) d^{3N}x},
$$

where \\(x\\) now represents all \\(3N\\) electron coordinates. To get this integral in
the standard Monte Carlo form, we define the quantity

$$
  E_L(x) = \frac{\hat{H}\psi(x)}{\psi(x)},
$$

called the *local energy* of \\(\psi\\), and the probability density

$$
  \rho(x) = \frac{|\psi(x)|^2}{\int |\psi(x')|^2 d^{3N}x}.
$$

The variational integral then becomes

$$
  \langle H \rangle = \langle E_L \rangle = \int \rho(x)E_L(x_) d^{3N}x.
$$

Now this is what we're looking for! All we now need to do is apply the
Metropolis-Hastings algorithm to sample electron configuraions \\(x_i\\), compute the
average of \\(E_L\\), \\(\bar{E}_L\\), on these points, and we have our integral! Additionally, we should
keep track of the variance of \\(E_L\\):

$$
  \mathrm{Var}[E_L] = \frac{1}{M} \sum_{i=1}^M (E_L(x_i) - \bar{E}_L)^2.
$$

Note that, if \\(\psi\\) is the exact ground state of \\(\hat{H}\\), we have
\\(E_L(x) = E_i\\) for all \\(x\\), and so the variance will vanish. Thus, the variance
of \\(E_L\\) is a direct measure of the quality of \\(\psi(x)\\).
