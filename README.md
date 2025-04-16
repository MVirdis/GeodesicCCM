Given a nonlinear smooth contraction metric M on a vector space and a differential controller dK, this tool allows to
  1. Find an approximate geodesic path for the metric M;
  2. integrate the differential controller along the path.

It uses the Clenshaw-Curtis Quadrature (CCQ) scheme to approximate integrals.

This is a MATLAB version of the accompanying code to the paper [Nonlinear Stabilization via Control Contraction Metrics: a
Pseudospectral Approach for Computing Geodesics](https://arxiv.org/pdf/1607.04340.pdf) by K. Leung and I. R. Manchester.
It does not implement the quasi-Newton BFGS method from scratch, but uses MATLAB's fmincon optimizer.
