# Analytical decay references and why each comparison is needed

The example's state is $(x,q)$: normalized concentration and exposure. Set $a=0.5$, $b=0.1$, $\kappa=0.1$, $\gamma=1$, $\beta=0.4$, and $\xi_p=3.2$. Let $t\in[t_0,t_f]=[2,12]$ minutes. These constants reproduce the revised Python verification objective while retaining a general nonzero initial time.

## State and integrated cost

For $\dot x=-x(a+bx)$, separation gives

$$
\frac{1}{a}[\log x-\log(a+bx)]=-\tau+C,\qquad\tau=t-t_0.
$$

Apply $x(t_0)=\xi$ and solve for $x$:

$$
r=e^{-a\tau},\qquad k=\frac ba(1-r),\qquad D=1+k\xi,
\qquad x=\frac{\xi r}{D}.
$$

Since $\dot D=b\xi r$, $\dot{\log D}/b=x$. Thus $q=\log D/b$ satisfies $\dot q=x$, $q(t_0)=0$. Integrating the concentration equation gives $x-\xi=-aq-b\int x^2dt$ and therefore

$$
C(t;\xi)=\frac\kappa2\int_{t_0}^t x^2\,ds
=\frac\kappa{2b}(\xi-x-aq).
$$

This analytical expression is a reference only. The C++ library computes the running cost by theta endpoint quadrature on its numerical physical trajectory, as required by the discrete formulation.

## Tangent and curvature

Differentiate at fixed elapsed time:

$$
x_\xi=\frac r{D^2},\quad q_\xi=\frac k{bD},\quad
x_{\xi\xi}=-\frac{2kr}{D^3},\quad q_{\xi\xi}=-\frac{k^2}{bD^2}.
$$

Then

$$
C_\xi=\frac\kappa{2b}(1-x_\xi-aq_\xi),
\qquad C_{\xi\xi}=\frac\kappa{2b}(-x_{\xi\xi}-aq_{\xi\xi}).
$$

Let $Q=q(t_f;3)$ and $r_q=q(t_f;\xi)-Q$. The verification objective has exact expressions

$$
J=C(t_f;\xi)+\frac\gamma2r_q^2+\frac\beta2(\xi-\xi_p)^2,
$$

$$
J'=C_\xi+\gamma r_q q_\xi+\beta(\xi-\xi_p),
\qquad J''=C_{\xi\xi}+\gamma(q_\xi^2+r_q q_{\xi\xi})+\beta.
$$

The independent reference optimizer brackets a zero of this analytical $J'$. The discrete Newton driver uses the library's adjoint and incremental adjoint, not these formulas. There is no claim of an elementary closed-form expression for the minimizing initial condition.

## Continuous costate

The continuous equations are

$$
\dot\lambda_x=(a+2bx)\lambda_x-\lambda_q-\kappa x,
\qquad\dot\lambda_q=0,
\qquad(\lambda_x(t_f),\lambda_q(t_f))=(0,\gamma r_q).
$$

An independent analytical expression follows by differentiating remaining cost. At time $t$, hold the present exposure fixed and regard $x(t)$ as a new initial concentration for the remaining duration $\Delta=t_f-t$. Write $S_x(\Delta;x)$ and $S_q(\Delta;x)$ for the sensitivity formulas above using that initial value and duration. Then

$$
\lambda_q=\gamma r_q,
\qquad\lambda_x(t)=\frac\kappa{2b}[1-S_x(\Delta;x(t))-aS_q(\Delta;x(t))]
+\gamma r_q S_q(\Delta;x(t)).
$$

This reference avoids numerically integrating a continuous adjoint. The library instead derives and solves the discrete adjoint; their difference measures temporal error after the documented endpoint convention is applied.

## Experimental design

At $\xi=2.4$, solve six meshes for each of $\theta=0.5$ and $0.65$. Fixed-input comparisons isolate temporal error in state, cost, tangent, adjoint, gradient, and Hessian. Reoptimize on every mesh and compare its initial condition separately to the analytical-objective minimizer. Mixing fitted trajectories and fixed-input reference trajectories would combine two different errors.

Finite differences test whether derivatives match the *discrete* maps. Analytical comparisons test whether those maps approximate the *continuous* problem. The global block assembly independently tests the coupling of forward, tangent, and adjoint steps. The two-control test exercises additional derivative paths absent from scalar decay. Together these provide stronger evidence than either an accurate forward plot or a gradient check alone.

For smooth stable trapezoidal integration, the local step defect is $O(h^3)$ and accumulated error is $O(h^2)$. Consistent cost quadrature and differentiated equations inherit this order in the stated smooth example. Near a nondegenerate isolated minimizer, $g_h(\xi^*)=O(h^2)$ and a Hessian bounded away from zero imply $\xi_h^*-\xi^*=O(h^2)$. Off-centering theta produces an $O(h^2)$ local defect and generally first-order global errors.

Finally, on a fixed grid compare Newton iterates with an independently bracketed zero of the *discrete* gradient. The local relation $e_{k+1}\approx0.1493e_k^2$ concerns optimizer iteration, not step size. Once the error reaches floating-point and nonlinear-solve tolerances, further iteration-order estimates are meaningless.
