# Appendix D. Linear algebra for inverse problems

## D.1 Symmetric positive definite matrices

A symmetric matrix \(\mathbf A\) is **positive definite** (SPD) if \(\mathbf x^\top\mathbf A\mathbf x>0\) for all \(\mathbf x\neq\mathbf 0\). Equivalent: all eigenvalues are
positive. SPD matrices are invertible, have a Cholesky factorisation \(\mathbf A = \mathbf L\mathbf L^\top\), and define the energy norm
\(\lVert\mathbf x\rVert_{\mathbf A} = \sqrt{\mathbf x^\top\mathbf A\mathbf x}\). The BoilFlux step matrix \(\mathbf A_n\) of (7.4) is SPD (Result 7.2).

**Gershgorin's theorem.** Every eigenvalue of a matrix lies in at least one disc centred at a diagonal entry \(A_{ii}\) with
radius \(\sum_{j\neq i}\lvert A_{ij}\rvert\). Used in Section 7.5 to bound the condition number.

## D.2 M-matrices and non-negative inverses

A matrix with positive diagonal and non-positive off-diagonal entries that is **strictly diagonally dominant**,
\(A_{ii}>\sum_{j\neq i}\lvert A_{ij}\rvert\), is a non-singular **M-matrix**.

!!! theorem "Lemma D.1"
    A non-singular M-matrix of the above kind has an entrywise non-negative inverse.

??? derivation "Proof"
    Write \(\mathbf A = \mathbf D(\mathbf I - \mathbf N)\) with \(\mathbf D = \operatorname{diag}(\mathbf A) > 0\) and \(\mathbf N = \mathbf I - \mathbf D^{-1}\mathbf A\). The matrix \(\mathbf N\) has zero
    diagonal and non-negative off-diagonal entries, and by strict diagonal dominance its row sums are below 1, so
    \(\lVert\mathbf N\rVert_\infty<1\). The Neumann series \((\mathbf I - \mathbf N)^{-1} = \sum_{k\ge0}\mathbf N^k\) converges and has non-negative terms. Hence
    \(\mathbf A^{-1} = (\mathbf I - \mathbf N)^{-1}\mathbf D^{-1}\ge0\). \(\square\)

This is the key step of the discrete maximum principle (Theorem 7.4).

## D.3 The singular value decomposition

Every real \(m\times n\) matrix has a factorisation

\[
\mathbf X = \mathbf U\mathbf S\mathbf V^\top = \sum_{i=1}^{r} s_i\,\mathbf u_i\mathbf v_i^\top ,
\tag{D.1}
\]

with orthonormal columns \(\mathbf u_i\in\mathbb R^m\), \(\mathbf v_i\in\mathbb R^n\), singular values \(s_1\ge\dots\ge s_r>0\) and \(r\) the rank. The
\(\mathbf v_i\) are eigenvectors of \(\mathbf X^\top\mathbf X\) with eigenvalues \(s_i^2\), and \(\mathbf u_i = \mathbf X\mathbf v_i/s_i\). Properties used in the book:

- \(\lVert\mathbf X\rVert_2 = s_1\); the condition number of a square invertible \(\mathbf X\) is \(\kappa = s_1/s_n\).
- The least-squares solution of minimum norm is \(\mathbf q = \sum_i(\mathbf u_i^\top\mathbf y/s_i)\,\mathbf v_i\) (8.7).
- Any function of \(\mathbf X^\top\mathbf X\) is diagonal in the basis \(\mathbf v_i\), which gives the filter factors of Chapter 9 and the
  posterior variances of Chapter 15.
- **Eckart–Young:** the best rank-\(k\) approximation of \(\mathbf X\) in the 2-norm is \(\sum_{i\le k}s_i\mathbf u_i\mathbf v_i^\top\), with error \(s_{k+1}\).

## D.4 The conjugate-gradient algorithm

For SPD \(\mathbf A\) and right-hand side \(\mathbf b\), starting from \(\mathbf x_0\), \(\mathbf r_0 = \mathbf b - \mathbf A\mathbf x_0\), \(\mathbf d_0 = \mathbf r_0\):

\[
\alpha_k = \frac{\mathbf r_k^\top\mathbf r_k}{\mathbf d_k^\top\mathbf A\mathbf d_k},\quad
\mathbf x_{k+1} = \mathbf x_k + \alpha_k\mathbf d_k,\quad
\mathbf r_{k+1} = \mathbf r_k - \alpha_k\mathbf A\mathbf d_k,\quad
\beta_{k+1} = \frac{\mathbf r_{k+1}^\top\mathbf r_{k+1}}{\mathbf r_k^\top\mathbf r_k},\quad
\mathbf d_{k+1} = \mathbf r_{k+1} + \beta_{k+1}\mathbf d_k .
\tag{D.2}
\]

Here \(\mathbf r_k = -\mathbf g_k\) is the negative gradient of \(\tfrac12\mathbf x^\top\mathbf A\mathbf x - \mathbf b^\top\mathbf x\). The iterate \(\mathbf x_k\) minimises the energy-norm error over
\(\mathbf x_0 + \mathcal K_k(\mathbf A, \mathbf r_0)\), which gives the convergence bound

\[
\lVert\mathbf x_k - \mathbf x^*\rVert_{\mathbf A} \le 2\left(\frac{\sqrt\kappa - 1}{\sqrt\kappa + 1}\right)^{k}\lVert\mathbf x_0 - \mathbf x^*\rVert_{\mathbf A}
\tag{D.3}
\]

through the Chebyshev polynomial of degree \(k\) on \([\lambda_\text{min}, \lambda_\text{max}]\). A **preconditioner** \(\mathbf P\approx\mathbf A\) replaces \(\kappa\)
by the condition number of \(\mathbf P^{-1}\mathbf A\); the Jacobi choice \(\mathbf P = \operatorname{diag}(\mathbf A)\) is used in BoilFlux.

**CGLS** applies (D.2) to the normal equations \(\mathbf X^\top\mathbf X\mathbf q = \mathbf X^\top\mathbf y\) without forming \(\mathbf X^\top\mathbf X\): each iteration needs one product
with \(\mathbf X\) and one with \(\mathbf X^\top\). In the quench problem these are one sensitivity solve and one adjoint solve (Chapter 10).
