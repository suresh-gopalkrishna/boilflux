# 12. Learning from data: networks, gradients and optimisation

<p class="bf-lead">Scientific machine learning uses the same mathematics as Part IV. Training a neural network
is a regularised least-squares problem, and back-propagation is the adjoint method of Chapter 10 applied to a
chain of layers. This chapter develops the foundations needed for physics-informed networks and neural
operators: approximation, gradients, optimisation, and the spectral bias that makes sharp boiling fronts hard
to learn.</p>

## 12.1 Learning as regularised approximation

Given data pairs \(\{(\mathbf x_i,\mathbf y_i)\}_{i=1}^N\) drawn from an unknown relation, choose a parametrised model
\(f_\theta\) and minimise the **empirical risk**

\[
\mathcal R_N(\theta) = \frac1N\sum_{i=1}^N\ell\big(f_\theta(\mathbf x_i),\mathbf y_i\big) + \frac\gamma2\lVert\theta\rVert^2 ,
\tag{12.1}
\]

usually with the squared loss \(\ell(\hat{\mathbf y},\mathbf y) = \lVert\hat{\mathbf y} - \mathbf y\rVert^2\). The second term (**weight decay**) is zeroth-order
Tikhonov regularisation on the parameters (Section 9.2). What we care about is the **expected risk** on new data;
the difference is the generalisation gap, estimated with held-out validation and test data. For physical
problems the test data should differ in the conditions that matter, for example an unseen heater temperature
(Chapter 14), and not only be a random subset of the training data.

## 12.2 Multilayer perceptrons and what they can represent

!!! definition "Definition 12.1 (multilayer perceptron)"
    A network with \(L\) layers maps \(\mathbf z_0 = \mathbf x\) to \(f_\theta(\mathbf x) = \mathbf W_L\mathbf z_{L-1} + \mathbf b_L\) through

    \[
    \mathbf z_l = \sigma\big(\mathbf W_l\,\mathbf z_{l-1} + \mathbf b_l\big),\qquad l = 1,\dots,L-1,
    \tag{12.2}
    \]

    with weights \(\mathbf W_l\), biases \(\mathbf b_l\) (together \(\theta\)) and an activation \(\sigma\) applied componentwise, for example
    \(\tanh\), the rectified linear unit \(\operatorname{ReLU}(s) = \max(0,s)\) or the Gaussian error linear unit (GELU).

!!! theorem "Theorem 12.1 (universal approximation; Leshno, Lin, Pinkus and Schocken, 1993)"
    Let \(\sigma\) be continuous and not a polynomial. Then for every continuous function \(g\) on a compact set
    \(K\subset\mathbb R^d\) and every \(\varepsilon>0\), there is a network with one hidden layer such that
    \(\sup_{\mathbf x\in K}\lvert f_\theta(\mathbf x) - g(\mathbf x)\rvert < \varepsilon\).

For ReLU in one dimension the proof is constructive. The piecewise-linear interpolant of \(g\) at nodes
\(x_0<x_1<\dots<x_m\) is

\[
g_m(x) = g(x_0) + \sum_{i=0}^{m-1} c_i\,\operatorname{ReLU}(x - x_i),\qquad
c_i = s_i - s_{i-1},\quad s_i = \frac{g(x_{i+1}) - g(x_i)}{x_{i+1}-x_i},\ s_{-1} = 0,
\tag{12.3}
\]

since each ReLU switches on a change of slope at its node. A uniformly continuous \(g\) is approximated
uniformly as the nodes are refined. The theorem guarantees existence only. It says nothing about how many
neurons are needed, or whether training finds them.

## 12.3 Back-propagation is the adjoint method

Write a network as the constraint chain \(\mathbf R_l = \mathbf z_l - \sigma(\mathbf a_l) = \mathbf 0\), \(\mathbf a_l = \mathbf W_l\mathbf z_{l-1} + \mathbf b_l\), with a loss
\(\ell(\mathbf z_L)\) (output layer linear for simplicity). This is exactly the setting of Section 10.2, with layers playing
the role of time steps. The Lagrangian is

\[
\mathcal L = \ell(\mathbf z_L) - \sum_{l=1}^{L}\boldsymbol\delta_l^\top\big(\mathbf z_l - \sigma(\mathbf a_l)\big) .
\tag{12.4}
\]

Setting \(\partial\mathcal L/\partial\mathbf z_l = \mathbf 0\): the state \(\mathbf z_l\) appears in constraint \(l\) and, through \(\mathbf a_{l+1}\), in constraint \(l+1\).

!!! theorem "Result 12.2 (back-propagation)"
    \[
    \boldsymbol\delta_L = \nabla\ell(\mathbf z_L),\qquad
    \boldsymbol\delta_l = \mathbf W_{l+1}^\top\big(\sigma'(\mathbf a_{l+1})\odot\boldsymbol\delta_{l+1}\big),\quad l = L-1,\dots,1,
    \tag{12.5}
    \]

    \[
    \frac{\partial\ell}{\partial\mathbf W_l} = \big(\sigma'(\mathbf a_l)\odot\boldsymbol\delta_l\big)\,\mathbf z_{l-1}^\top,\qquad
    \frac{\partial\ell}{\partial\mathbf b_l} = \sigma'(\mathbf a_l)\odot\boldsymbol\delta_l .
    \tag{12.6}
    \]

    (For a linear output layer set \(\sigma'(\mathbf a_L) = \mathbf 1\).)

Compare with (10.12): a backward recursion, started from the derivative of the loss at the last stage, that
applies the **transposed** Jacobian of each stage, followed by a gradient formula that pairs the adjoint with the
input of the stage. The memory requirement is the same as well: the forward states must be stored. Automatic
differentiation in PyTorch or JAX implements (12.5)–(12.6) for any program. It costs a small multiple of one
forward evaluation, however many parameters there are.

**Forward versus reverse mode.** For a map \(\mathbb R^n\to\mathbb R^m\), forward-mode differentiation costs one pass per
input direction, reverse mode one pass per output. Training (one scalar loss, millions of parameters) needs reverse
mode. PINNs also need derivatives of the network output with respect to its few input coordinates (Chapter 13),
for which nested reverse mode or forward mode is used.

## 12.4 Optimisation

**Gradient descent** \(\theta_{k+1} = \theta_k - \eta\nabla\mathcal R(\theta_k)\) decreases a function whose gradient is \(\Lambda\)-Lipschitz if
\(\eta<2/\Lambda\): by Taylor's theorem, \(\mathcal R(\theta_{k+1})\le\mathcal R(\theta_k) - \eta(1 - \eta\Lambda/2)\lVert\nabla\mathcal R\rVert^2\). Its speed on an
ill-conditioned problem is limited exactly as in Section 7.5, which is why inputs and outputs are normalised.
**Stochastic** gradient descent estimates \(\nabla\mathcal R\) from a mini-batch, which is cheaper and adds noise that helps to escape
poor stationary points.

**Adam** keeps running averages of the gradient \(\mathbf g_k\) and of its square,

\[
\mathbf m_k = \beta_1\mathbf m_{k-1} + (1-\beta_1)\mathbf g_k,\qquad
\mathbf v_k = \beta_2\mathbf v_{k-1} + (1-\beta_2)\,\mathbf g_k\odot\mathbf g_k,\qquad \mathbf m_0 = \mathbf v_0 = \mathbf 0 ,
\tag{12.7}
\]

and updates \(\theta_{k+1} = \theta_k - \eta\,\hat{\mathbf m}_k/(\sqrt{\hat{\mathbf v}_k}+\epsilon)\). The **bias correction**
\(\hat{\mathbf m}_k = \mathbf m_k/(1-\beta_1^k)\) follows from unrolling the recursion: \(\mathbf m_k = (1-\beta_1)\sum_{j=1}^k\beta_1^{k-j}\mathbf g_j\), so for
stationary gradients \(\mathbb E\,\mathbf m_k = (1-\beta_1^k)\,\mathbb E\,\mathbf g\). Without correction the first steps would be too small. The
division by \(\sqrt{\hat{\mathbf v}}\) gives every parameter a step of roughly \(\eta\), a diagonal preconditioner. **AdamW** applies weight
decay directly to the parameters instead of through the gradient, so that it is not rescaled by the
preconditioner. The neural operator of Chapter 14 is trained with AdamW, a one-cycle learning-rate schedule and
gradient clipping.

## 12.5 Spectral bias

Training a wide network behaves, to first order, like a linear model in its parameters. Linearise about the initial
parameters, \(f_\theta(\mathbf x)\approx f_{\theta_0}(\mathbf x) + \nabla_\theta f_{\theta_0}(\mathbf x)^\top(\theta - \theta_0)\), and run gradient flow
\(\dot\theta = -\eta\nabla_\theta\mathcal R\) on the squared loss. The vector of residuals at the training points, \(\mathbf r(t)\), then obeys

\[
\dot{\mathbf r} = -\eta\,\boldsymbol\Theta\,\mathbf r,\qquad
\Theta_{ij} = \nabla_\theta f_{\theta_0}(\mathbf x_i)^\top\nabla_\theta f_{\theta_0}(\mathbf x_j),
\qquad\Longrightarrow\qquad
\mathbf r(t) = \sum_k e^{-\eta\lambda_k t}\,(\boldsymbol\phi_k^\top\mathbf r_0)\,\boldsymbol\phi_k ,
\tag{12.8}
\]

with the eigenpairs \((\lambda_k,\boldsymbol\phi_k)\) of the **neural tangent kernel** \(\boldsymbol\Theta\). For standard networks and activations, the
eigenvalues decrease rapidly with the frequency content of the eigenvectors. Smooth components of the target are
learned first and sharp ones last, exactly like the iterative regularisation of Section 9.4.

This **spectral bias** is a form of regularisation, but it is a problem when the target is genuinely sharp: wetting
fronts, vapour–liquid interfaces and the rewetting surge are steep in space and time. Remedies include mapping the
inputs through **Fourier features** \(\gamma(\mathbf x) = [\cos(2\pi\mathbf B\mathbf x),\ \sin(2\pi\mathbf B\mathbf x)]\) with a random matrix \(\mathbf B\), which spreads
the kernel spectrum, and architectures that work in a spectral basis directly (Chapter 14).

## Exercises

!!! exercise "12.1"
    Verify (12.3) for \(g(x) = x^2\) on \([0,1]\) with three equal intervals, and bound the uniform error.

!!! exercise "12.2"
    Derive (12.5)–(12.6) from (12.4) by differentiating with respect to \(\mathbf z_l\), \(\mathbf W_l\) and \(\mathbf b_l\).

!!! exercise "12.3"
    Show that for a linear model \(f_\theta(\mathbf x) = \theta^\top\boldsymbol\varphi(\mathbf x)\) the kernel \(\boldsymbol\Theta\) does not change during training, so (12.8)
    is exact, and relate its eigen-decomposition to the SVD filter factors of Landweber iteration (9.10).
