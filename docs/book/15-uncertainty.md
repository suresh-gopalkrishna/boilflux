# 15. Bayesian inverse problems and uncertainty

<p class="bf-lead">A reconstruction without an error bar cannot be compared with a model or with another
experiment. The Bayesian view turns regularisation into prior information and gives the uncertainty that the
data leave. This chapter derives the Gaussian posterior of the linear problem, shows that its mean is the
Tikhonov solution, explains what the posterior variance and the resolution matrix say about a quench
reconstruction, and describes how to estimate uncertainty for the nonlinear problem and for learned models.</p>

## 15.1 Bayes' theorem for inverse problems

Treat the unknown flux \(\mathbf q\) as a random vector. Before the experiment, a **prior** density \(\pi_\text{pr}(\mathbf q)\) states what
is plausible. The measurement model gives the **likelihood** \(\pi(\mathbf y\mid\mathbf q)\). Bayes' theorem combines them into the
**posterior**:

\[
\pi_\text{post}(\mathbf q\mid\mathbf y) = \frac{\pi(\mathbf y\mid\mathbf q)\,\pi_\text{pr}(\mathbf q)}{\pi(\mathbf y)} \ \propto\ \pi(\mathbf y\mid\mathbf q)\,\pi_\text{pr}(\mathbf q) .
\tag{15.1}
\]

The posterior is the complete answer to the inverse problem: every estimate and error bar is a summary of it.

## 15.2 The linear Gaussian case

Let \(\mathbf y = \mathbf X\mathbf q + \mathbf e\) with Gaussian noise \(\mathbf e\sim\mathcal N(\mathbf 0, \sigma^2\mathbf I)\) and a Gaussian prior
\(\mathbf q\sim\mathcal N(\mathbf q_0, \boldsymbol\Gamma_\text{pr})\). Then

\[
-2\ln\pi_\text{post} = \frac{\lVert\mathbf X\mathbf q - \mathbf y\rVert^2}{\sigma^2} + (\mathbf q - \mathbf q_0)^\top\boldsymbol\Gamma_\text{pr}^{-1}(\mathbf q - \mathbf q_0) + \text{const} .
\]

The right side is quadratic in \(\mathbf q\). Collecting the quadratic and linear terms,

\[
-2\ln\pi_\text{post} = \mathbf q^\top\mathbf P\mathbf q - 2\mathbf q^\top\big(\mathbf X^\top\mathbf y/\sigma^2 + \boldsymbol\Gamma_\text{pr}^{-1}\mathbf q_0\big) + \text{const},
\qquad \mathbf P = \frac{\mathbf X^\top\mathbf X}{\sigma^2} + \boldsymbol\Gamma_\text{pr}^{-1} ,
\]

and completing the square identifies a Gaussian:

!!! theorem "Result 15.1 (Gaussian posterior)"
    \[
    \pi_\text{post} = \mathcal N(\mathbf q_\text{post},\,\boldsymbol\Gamma_\text{post}),\qquad
    \boldsymbol\Gamma_\text{post} = \mathbf P^{-1},\qquad
    \mathbf q_\text{post} = \mathbf P^{-1}\Big(\frac{\mathbf X^\top\mathbf y}{\sigma^2} + \boldsymbol\Gamma_\text{pr}^{-1}\mathbf q_0\Big).
    \tag{15.2}
    \]

    For the prior \(\mathbf q_0 = \mathbf 0\), \(\boldsymbol\Gamma_\text{pr} = (\sigma^2/\gamma)\,\mathbf I\), the posterior mean (and maximum) is the Tikhonov solution (9.3),
    and \(\boldsymbol\Gamma_\text{post} = \sigma^2(\mathbf X^\top\mathbf X + \gamma\mathbf I)^{-1}\).

Tikhonov regularisation is therefore the **maximum a posteriori** (MAP) estimate under a Gaussian prior of standard
deviation \(\sigma/\sqrt\gamma\) per flux value. Higher-order Tikhonov corresponds to a prior on derivatives (a smoothness prior).
The regularisation parameter acquires a meaning: it is the ratio of noise variance to prior variance.

## 15.3 What the posterior says

In the SVD basis, \(\mathbf X^\top\mathbf X + \gamma\mathbf I = \sum_i(s_i^2 + \gamma)\mathbf v_i\mathbf v_i^\top\), so the posterior variance along the singular vector
\(\mathbf v_i\) is

\[
\operatorname{Var}\big(\mathbf v_i^\top\mathbf q\mid\mathbf y\big) = \frac{\sigma^2}{s_i^2 + \gamma} .
\tag{15.3}
\]

Components with \(s_i^2\gg\gamma\) are determined by the data, with variance \(\sigma^2/s_i^2\). Components with \(s_i^2\ll\gamma\) are
not: their variance is the prior variance \(\sigma^2/\gamma\). The posterior is honest about the ill-posedness of Chapter 8: the
uncertainty of fast and fine-scale flux content equals what we assumed before measuring.

**Resolution.** The posterior mean is a filtered version of the truth plus filtered noise:

\[
\mathbf q_\text{post} = \mathbf R\,\mathbf q^\dagger + (\mathbf X^\top\mathbf X + \gamma\mathbf I)^{-1}\mathbf X^\top\mathbf e,
\qquad
\mathbf R = (\mathbf X^\top\mathbf X + \gamma\mathbf I)^{-1}\mathbf X^\top\mathbf X = \sum_i f_i\,\mathbf v_i\mathbf v_i^\top .
\tag{15.4}
\]

Row \(j\) of the **resolution matrix** \(\mathbf R\) is the **averaging kernel** of flux value \(j\): the weights with which the true flux
at all points and times contributes to the estimate there. Its width is the actual space and time resolution of the
reconstruction, and it is the right quantity to quote next to a peak value such as \(q_\text{max}\).

## 15.4 The nonlinear problem: the Laplace approximation

For temperature-dependent properties the posterior is not Gaussian. The **Laplace approximation** replaces it by a
Gaussian centred at the MAP estimate \(\mathbf q^*\), with covariance equal to the inverse Hessian of \(-\ln\pi_\text{post}\) there. In the
Gauss–Newton approximation,

\[
\boldsymbol\Gamma_\text{post}^{-1}\approx\frac{\mathbf X(\mathbf q^*)^\top\mathbf X(\mathbf q^*)}{\sigma^2} + \boldsymbol\Gamma_\text{pr}^{-1},
\tag{15.5}
\]

with \(\mathbf X(\mathbf q^*)\) the sensitivity matrix at the solution. It is never formed. A product
\(\mathbf X^\top\mathbf X\mathbf v\) costs one sensitivity solve (10.14) for \(\mathbf X\mathbf v\) and one adjoint solve (10.12) for \(\mathbf X^\top(\cdot)\). A Lanczos
iteration with these products finds the dominant eigenpairs of the data term; since only a few dozen
directions are informed by the data (Figure 8.1), a low-rank update of the prior covariance represents
\(\boldsymbol\Gamma_\text{post}\) accurately at the cost of a few dozen pairs of solves.

## 15.5 Sampling by perturbing the data

A simple and general way to explore the posterior is to repeat the inversion with perturbed inputs.

!!! theorem "Result 15.2 (randomise-then-optimise, linear Gaussian case)"
    Draw \(\mathbf e_j\sim\mathcal N(\mathbf 0,\sigma^2\mathbf I)\) and \(\boldsymbol\eta_j\sim\mathcal N(\mathbf 0,\boldsymbol\Gamma_\text{pr})\), and solve

    \[
    \mathbf q_j = \arg\min_{\mathbf q}\ \frac{\lVert\mathbf X\mathbf q - (\mathbf y + \mathbf e_j)\rVert^2}{\sigma^2} + (\mathbf q - \mathbf q_0 - \boldsymbol\eta_j)^\top\boldsymbol\Gamma_\text{pr}^{-1}(\mathbf q - \mathbf q_0 - \boldsymbol\eta_j).
    \tag{15.6}
    \]

    Then the \(\mathbf q_j\) are exact, independent samples of the posterior (15.2).

??? derivation "Proof"
    The minimiser is \(\mathbf q_j = \mathbf P^{-1}\big(\mathbf X^\top(\mathbf y+\mathbf e_j)/\sigma^2 + \boldsymbol\Gamma_\text{pr}^{-1}(\mathbf q_0+\boldsymbol\eta_j)\big)\), an affine function of Gaussian
    variables, hence Gaussian. Its mean is \(\mathbf q_\text{post}\). The two random terms are independent with covariances
    \(\mathbf X^\top\mathbf X/\sigma^2\) and \(\boldsymbol\Gamma_\text{pr}^{-1}\), which add up to \(\mathbf P\). So \(\operatorname{Cov}(\mathbf q_j) = \mathbf P^{-1}\mathbf P\mathbf P^{-1} = \mathbf P^{-1}\). \(\square\)

For the nonlinear problem the samples are approximate, but the procedure is the same: rerun
`AdjointCGM` on data with fresh noise draws. With iterative regularisation (no explicit prior), perturbing only the
data and keeping the discrepancy stop gives a practical estimate of the **noise-induced** spread. It does not capture
the bias of regularisation, which the averaging kernels of Section 15.3 describe.

## 15.6 Model error

Noise is often not the dominant error. The model itself carries uncertainty:

| Source | Effect on the reconstructed flux | How to account for it |
|---|---|---|
| thermophysical properties \(k, \rho c_p\) | flux scales roughly with \(\rho c_p\) (energy) and with \(k\) (through-thickness gradient) | propagate property uncertainty through repeated inversions |
| IR calibration and emissivity | temperature bias passes into the flux through the energy balance | calibrate in situ; include as a systematic error |
| initial temperature field | early-time flux biased | measure \(T_0\) on the IR face before the quench |
| unmodelled latent heat or contact resistance | spurious flux features at transformation temperatures | effective \(c_p\), Section 1.5 |
| discretisation | small by Chapter 6–7 estimates | grid refinement study |

The **Bayesian approximation-error** approach estimates the statistics of the model error by comparing an accurate and
the inversion model on prior samples and adds them to the noise covariance. BF-Bench's refined-grid data generation
(Section 5.5) is a first step in this direction.

## 15.7 Uncertainty of learned models

For neural networks the posterior over weights is intractable. Practical approximations are **deep ensembles**
(train several networks from different initialisations and data orders; the spread estimates the epistemic uncertainty),
Monte-Carlo dropout, and variational or Laplace approximations of Bayesian neural networks. Whatever the method,
check **calibration**: the fraction of held-out values inside a nominal 90 % interval should be close to 90 %. A surrogate
used inside an inverse method (Chapter 16) must report when it is outside its training distribution, because there
its error bars are least reliable.

## 15.8 What to report with a reconstruction

1. The noise level used for stopping and how it was obtained (given or estimated by (11.5)).
2. The reliable time window (Section 10.7) and the observed region (mask).
3. The spatial and temporal resolution: widths of averaging kernels, or the plate thickness and \(t_d\) as a rule of thumb.
4. Noise-induced spread from repeated inversions with perturbed data.
5. Sensitivity to the property tables and the emissivity.

## Exercises

!!! exercise "15.1"
    Show that the posterior mean (15.2) can also be written
    \(\mathbf q_\text{post} = \mathbf q_0 + \boldsymbol\Gamma_\text{pr}\mathbf X^\top(\mathbf X\boldsymbol\Gamma_\text{pr}\mathbf X^\top + \sigma^2\mathbf I)^{-1}(\mathbf y - \mathbf X\mathbf q_0)\), and explain when this form is cheaper.

!!! exercise "15.2"
    For the 1D example of Section 8.4 with \(\gamma\) from the discrepancy principle, compute the posterior standard
    deviation of the flux at \(t = 0.3\) s from (15.3) and compare with the reconstruction error in Figure 9.2.

!!! exercise "15.3"
    Show that the trace of \(\mathbf R\) in (15.4) equals \(\sum_i f_i\), and interpret it as the number of flux parameters that the
    data determine.
