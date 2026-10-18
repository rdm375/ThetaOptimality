# The theta-consistent discrete optimality system

We begin with the forward update rather than the residual notation used internally by Newton. Let $\mathbf x(t)\in\mathbb R^N$, $\boldsymbol\xi\in\mathbb R^P$, $t\in[t_0,t_f]$, and

$$
\dot{\mathbf x}=\mathbf f(t,\mathbf x,\boldsymbol\xi),
\qquad \mathbf x(t_0)=\mathbf X_0(\boldsymbol\xi).
$$

Both the initial map and dynamics can depend on controls. Time nodes and theta weights are prescribed and held fixed during differentiation.

## Forward update and objective

$$
\mathbf x_{n+1}=\mathbf x_n+h_n[(1-\theta_n)\mathbf f_n+\theta_n\mathbf f_{n+1}],
\qquad \mathbf x_0=\mathbf X_0(\boldsymbol\xi).
$$

$$
J_h=\sum_{n=0}^{K-1}h_n[(1-\theta_n)c_n+\theta_n c_{n+1}]
+\Phi(t_f,\mathbf x_K,\boldsymbol\xi)+\mathcal R(\boldsymbol\xi).
$$

The endpoint weights combine into $w_0=h_0(1-\theta_0)$, $w_K=h_{K-1}\theta_{K-1}$, and

$$
w_n=h_{n-1}\theta_{n-1}+h_n(1-\theta_n),\quad 0<n<K.
$$

Thus $J_h=\sum_n w_n c_n+\Phi+\mathcal R$. This matters on nonuniform grids: an interior weight is generally not just one interval length. With isolated observations, add explicit costs $\sum_j\ell_j(t_{n_j},\mathbf x_{n_j},\boldsymbol\xi)$, without multiplying them by $w_{n_j}$. Their state gradients add to $\mathbf q_{n_j}$, their direct control gradients add to the reduced gradient, and their directional gradient derivatives add to the incremental adjoint. The same derivation otherwise applies. Pure sampled-data assimilation has zero running cost; theta-consistent quadrature remains available for genuine running terms and is tested in the separate running-cost example.

For implementation, write the same update as $\mathbf R_n=\mathbf x_{n+1}-\mathbf x_n-h_n[(1-\theta_n)\mathbf f_n+\theta_n\mathbf f_{n+1}]=0$. Its derivative blocks are $\partial_{\mathbf x_{n+1}}\mathbf R_n=\mathbf B_n$ and $\partial_{\mathbf x_n}\mathbf R_n=-\mathbf C_n$, with

$$
\mathbf B_n=\mathbf I-h_n\theta_n\nabla_{\mathbf x}\mathbf f_{n+1},
\qquad \mathbf C_n=\mathbf I+h_n(1-\theta_n)\nabla_{\mathbf x}\mathbf f_n.
$$

The library uses damped Newton to solve this algebraic equation; the residual is an implementation form of the theta update, not a different physical model.

## Tangent linear model

For a control direction $\mathbf v$, let $\mathbf s_n=D_{\boldsymbol\xi}\mathbf x_n[\mathbf v]$. Define

$$
\mathbf K_n=h_n[(1-\theta_n)\nabla_{\boldsymbol\xi}\mathbf f_n+\theta_n\nabla_{\boldsymbol\xi}\mathbf f_{n+1}].
$$

Differentiating the discrete update gives

$$
\mathbf s_0=\mathbf G\mathbf v,\quad \mathbf G=\nabla_{\boldsymbol\xi}\mathbf X_0,
\qquad \mathbf B_n\mathbf s_{n+1}=\mathbf C_n\mathbf s_n+\mathbf K_n\mathbf v.
$$

This derivative is taken after discretization. There is no independently time-discretized continuous sensitivity equation in the primary implementation.

## Discrete adjoint and reduced gradient

Let

$$
\mathbf q_n=w_n\nabla_{\mathbf x}c_n+\mathbf1_{n=K}\nabla_{\mathbf x}\Phi,
\qquad \boldsymbol\lambda_K=\mathbf q_K.
$$

Sweep backward:

$$
\mathbf B_n^T\mathbf z_n=\boldsymbol\lambda_{n+1},
\qquad\boldsymbol\lambda_n=\mathbf q_n+\mathbf C_n^T\mathbf z_n.
$$

The reduced gradient is

$$
\nabla J_h=\nabla\mathcal R
+\sum_{n=0}^K w_n\nabla_{\boldsymbol\xi}c_n
+\nabla_{\boldsymbol\xi}\Phi
+\mathbf G^T\boldsymbol\lambda_0
+\sum_{n=0}^{K-1}\mathbf K_n^T\mathbf z_n.
$$

All cost control derivatives here hold the state fixed; state dependence has already been handled by the adjoint. Including it again would double-count it. When controls enter only as $\mathbf x_0=\boldsymbol\xi$, the formula reduces to $\nabla J_h=\boldsymbol\lambda_0+\nabla\mathcal R$.

The scalar product of this gradient with $\mathbf v$ must agree with direct tangent differentiation:

$$
D J_h[\mathbf v]=\nabla\mathcal R^T\mathbf v
+\sum_n w_n[\nabla_{\mathbf x}c_n^T\mathbf s_n+\nabla_{\boldsymbol\xi}c_n^T\mathbf v]
+\nabla_{\mathbf x}\Phi^T\mathbf s_K+\nabla_{\boldsymbol\xi}\Phi^T\mathbf v.
$$

## Why this is an optimality system

Stack the initial constraint and step residuals into $\mathbf F(\mathbf X,\boldsymbol\xi)=0$. With Lagrangian $\mathcal L=J_h-\mathbf p^T\mathbf F$, the first-order conditions are primal feasibility, $\nabla_{\mathbf X}\mathcal L=0$, and $\nabla_{\boldsymbol\xi}\mathcal L=0$. These are the discrete optimality or KKT equations for an interior optimum.

The tangent equation is linearized feasibility. It accompanies the KKT system and is used to differentiate it; it is not a separate stationarity condition. The library exploits the forward/backward block structure rather than assembling and factoring a full nonlinear KKT matrix. Tests independently assemble a small global system to verify this equivalence.

## Incremental adjoint and exact Hessian-vector product

Differentiate the same adjoint along $(\mathbf s_n,\mathbf v)$. All quantities below are directional derivatives of the discrete objects already defined:

$$
\delta\boldsymbol\lambda_K=\delta\mathbf q_K,
$$

$$
\mathbf B_n^T\delta\mathbf z_n
=\delta\boldsymbol\lambda_{n+1}-(\delta\mathbf B_n)^T\mathbf z_n,
$$

$$
\delta\boldsymbol\lambda_n=\delta\mathbf q_n
+(\delta\mathbf C_n)^T\mathbf z_n+\mathbf C_n^T\delta\mathbf z_n.
$$

Here $\delta\mathbf B_n=-h_n\theta_n\delta(\nabla_{\mathbf x}\mathbf f_{n+1})$ and $\delta\mathbf C_n=h_n(1-\theta_n)\delta(\nabla_{\mathbf x}\mathbf f_n)$. Differentiating the reduced gradient yields

$$
\begin{aligned}
\nabla^2J_h\,\mathbf v={}&\nabla^2\mathcal R\,\mathbf v
+\sum_n w_n\delta(\nabla_{\boldsymbol\xi}c_n)
+\delta(\nabla_{\boldsymbol\xi}\Phi)\\
&+(\delta\mathbf G)^T\boldsymbol\lambda_0+\mathbf G^T\delta\boldsymbol\lambda_0\\
&+\sum_n[(\delta\mathbf K_n)^T\mathbf z_n+\mathbf K_n^T\delta\mathbf z_n].
\end{aligned}
$$

This includes initial-map curvature, mixed state/control terms, dynamics curvature, and cost curvature. No term is dropped as in a Gauss–Newton approximation. Directional derivative callbacks supply these contractions without storing full third-order tensors.

## Nodal adjoint convention and continuous comparison

The raw quadrature-source nodal adjoint includes $h_{n-1}\theta_{n-1}\nabla_{\mathbf x}c_n$ from the preceding interval. A remaining-cost derivative at node $n$ excludes this already accumulated contribution:

$$
\boldsymbol\lambda_n^{\mathrm{future}}=\boldsymbol\lambda_n
-h_{n-1}\theta_{n-1}\nabla_{\mathbf x}c_n,\quad n\ge1,
\qquad\boldsymbol\lambda_0^{\mathrm{future}}=\boldsymbol\lambda_0.
$$

That is the convention used to compare against the continuous costate. The raw nodal adjoint remains the primary output and the gradient uses it. Continuous costates, nodal adjoints, and step multipliers should never be interchanged solely because all are called adjoints.

With smooth data and stable propagation, constant $\theta=1/2$ gives second-order state and cost convergence. Consistent differentiation then gives the corresponding sensitivity and properly interpreted adjoint convergence in this smooth example. The decay experiment measures those orders explicitly, including the reduced Hessian and minimizer. These are analytical comparisons, separate from finite-difference checks that establish exact differentiation of the discrete problem.
