# Distributions of the extent and mass of high excursion components of Gaussian random fields: verification of asymptotic distributions by simulation

R code for verification of theoretical distributions in ''Distributions of the extent and mass of high excursion components of Gaussian random fields'' by Dan Cheng, Zhibing He and Haoran Shi. 

The paper is still in progress and its link can not be provided for the moment.

## Main Body

Let $X(t),t\in \mathbb{R}^N$ be a smooth Gaussian random field, and let $A_u(t_0),u>0$ be the component of the excursion set $\{t\in \mathbb{R}^N:X(t)\geq u\}$ containing a peak at $t_0$. We will verify the asymptotic distributions for the extent (size of $A_u(t_0)$) and mass (volume of $A_u(t_0)$ under $X(t)$). In the paper, we consider the following two cases:

- $X(t)$ is isotropic, the random matric theory, such as GOI and GOE, can be applied 
- $X(t)$ is not isotropic but is stationary, Hessian of the field conditioned on a large local maximum will be derived

### 1. The 1D case

For convenience, we let $X(t)$ be a stationary centered unit-variance smooth Gaussian process, and define
$$
\lambda_1 = {\rm Var}(X'(0)), \quad \lambda_2 = {\rm Var}(X''(0)).
$$
Let $\phi(x)$ and $\Phi(x)$ be the pdf and cdf of the standard normal distribution respectively, and we define $\Psi(x) = 1 - \Phi(x)$.

### 1.1 Distribution of excursion length

The length of a excursion component above threshold $u$ is given by
$$
\textnormal{excursion length} = 2\sqrt{\frac{2(X(0)-u)}{-X''(0)}}
$$
The tail probability that the excursion length is greater than $\ell$ is given by
$$
\begin{align}\tag{1}
F_u(l) &= \mathbb{P}[\textnormal{excursion length}>\ell | X(0)>u,X'(0)=0,X''(0)<0] \\
       &= \frac{\frac{\sqrt{\lambda_2}}{\sqrt{2\pi}}\Psi\left( \frac{\sqrt{\lambda_2}}{\sqrt{\lambda_2-\lambda_1^2}}u\right) - \frac{\lambda_2\ell^2-8\lambda_1}{\sqrt{\lambda_2\ell^4-16\lambda_1\ell^2+64}}\phi\left(\frac{8u}{\sqrt{\lambda_2\ell^4-16\lambda_1\ell^2+64}}\right)\Psi\left( \frac{(\lambda_2\ell^2-8\lambda_1)u}{\sqrt{\lambda_2-\lambda_1^2}\sqrt{\lambda_2\ell^4-16\lambda_1\ell^2+64}}\right)}{\frac{\sqrt{\lambda_2}}{\sqrt{2\pi}}\Psi\left( \frac{\sqrt{\lambda_2}}{\sqrt{\lambda_2-\lambda_1^2}}u\right) + \lambda_1\phi(u)\Psi\left( \frac{-\lambda_1u}{\sqrt{\lambda_2-\lambda_1^2}}\right)}
\end{align}
$$

#### 1.2 Distribution of excursion area

The area of a excursion component above threshold $u$ is given by
$$
\textnormal{excursion area} = \frac{4\sqrt{2}}{3}\sqrt{\frac{(X(0)-u)^3}{-X''(0)}}
$$
The tail probability that the excursion area is greater than $V$ is given by
$$
\begin{align}\tag{2}
G_u(V) &= \mathbb{P}[\textnormal{excursion area} > V | X(0)>u,X'(0)=0,X''(0)<0] \\
       &= \left[\frac{\sqrt{\lambda_2}}{\sqrt{2\pi}}\Psi\left( \frac{\sqrt{\lambda_2}}{\sqrt{\lambda_2-\lambda_1^2}}u\right) + \lambda_1\phi(u)\Psi\left( \frac{-\lambda_1u}{\sqrt{\lambda_2-\lambda_1^2}}\right)\right]^{-1}\\
			&\quad \times \frac{2048}{81\sqrt{2\pi}\sqrt{\lambda_2-\lambda_1^2}}\int_V^\infty\frac{1}{y^5}\int_0^\infty  x^6 \exp\{-\frac{\left[\lambda_1(x+u)-\frac{32x^3}{9y^2}\right]^2}{2(\lambda_2-\lambda_1^2)}\}\phi(x+u)dx\,dy
\end{align}
$$

### 2. The 2D case

#### 2.1 Distribution of excursion extent

The extent of a excursion component above threshold $u$ is given by
$$
\textnormal{exursion extent} = \pi\frac{2(X(0)-u)}{\sqrt{{\rm det}\nabla^2 X(0)}}.
$$
The tail probability that the excursion extent greater than $\ell$ is given by
$$
\begin{align}\tag{3}
F_u(\ell) &= \mathbb{P}\left[\textnormal{exursion extent}>\ell \Big| X(0)>u, \nabla X(0)=0, \nabla^2 X(0)\prec 0 \right] \\
& = \frac{\int_u^\infty\mathbb{E} \left[(\Lambda_2-\Lambda_1)(\Lambda_1-\frac{\kappa x}{\sqrt{2}})(\Lambda_2-\frac{\kappa x}{\sqrt{2}})\mathbb{1}_1 \mathbb{1}_2\right] \phi(x)dx}{\int_u^\infty \mathbb{E} \left[(\Lambda_2-\Lambda_1)(\Lambda_1-\frac{\kappa x}{\sqrt{2}})(\Lambda_2 - \frac{\kappa x}{\sqrt{2}})\mathbb{1}_1\right]\phi(x)dx},
\end{align}
$$
where $\Lambda_1$ and $\Lambda_2$ are two i.i.d. standard normal variables, $\kappa = 1$, and 
$$
\begin{align}
\mathbb{1}_2 &= \mathbb{1}\{\frac{2\pi(x-u)}{\sqrt{(\Lambda_1-\kappa x/\sqrt{2})(\Lambda_2-\kappa x/\sqrt{2})/(2\gamma^4)}}>\ell\}, \\
\mathbb{1}_1 &=  \mathbb{1}\{\Lambda_1<\Lambda_2<\kappa x/\sqrt{2}\}.
\end{align}
$$


#### 2.2 Distribution of excursion volume

The volume of a excursion component above threshold $u$ is given by
$$
\textnormal{excursion volume} = \pi\frac{(X(0)-u)^2}{\sqrt{{\rm det}\nabla^2 X(0)}}.
$$
The tail probability that the excursion volume is greater than $v$ is given by 
$$
\begin{align}\tag{4}
F_u(v) &= \mathbb{P}\left[\textnormal{excursion volume} > v \Big| X(0)>u, \nabla X(0)=0, \nabla^2 X(0)\prec 0 \right], \\
&= \frac{\int_u^\infty\mathbb{E} \left[(\Lambda_2-\Lambda_1)(\Lambda_1-\frac{\kappa x}{\sqrt{2}})(\Lambda_2-\frac{\kappa x}{\sqrt{2}})\mathbb{1}_1 \mathbb{1}_3\right]\phi(x)dx}{\int_u^\infty \mathbb{E} \left[(\Lambda_2-\Lambda_1)(\Lambda_1-\kappa x/\sqrt{2})(\Lambda_2-\kappa x/\sqrt{2})\mathbb{1}_1\right]\phi(x)dx},
\end{align}
$$
where 
$$
\mathbb{1}_3 = \mathbb{1}\{\frac{\pi(x-u)^2}{\sqrt{(\Lambda_1-\kappa x/\sqrt{2})(\Lambda_2-\kappa x/\sqrt{2})/(2\gamma^4)}}>v\}.
$$

### 3. Joint Distribution

#### 3.1 Joint distribution of peak height and excursion length

#### 3.2 Joint distribution of peak height and excursion length



















