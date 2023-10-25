# DryTooling.FiniteVolumes theory

```@contents
Pages = ["theory.md"]
Depth = 3
```

## Heat conduction

### Cylindrical coordinates 1-D

Heat equation formulated with temperature as dependent variable is stated as:

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=\nabla\cdotp{}(k\nabla{}T)
```

For computing the heating dynamics in a cylinder, using the definition of divergence in cylindrical coordinates and using the gradient expansion over the radius we have

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=
\frac{1}{r}\frac{\partial}{\partial{}r}
\left(rk\frac{\partial{}T}{\partial{}r}\right)
```

To proceed with the finite volume discretization we perform the integration of both sides of the equation over the relevant variables. The order of integration is chosen according to the nature of the derivative term, as discussed by Patankar [Patankar1980](@cite). Care must be taken in the definition of the space integration, which is non-trivial in cylindrical coordinates systems and must be carried over the differential volume ``dV``.

```math
\int_{V}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}dtdV=
\int_{0}^{\tau}\int_{V}
\frac{1}{r}\frac{\partial}{\partial{}r}
\left(rk\frac{\partial{}T}{\partial{}r}\right)dVdt
```

This differential volume is given by ``dV=rdr{}d\theta{}dz``. Since the problem is specified to be symmetric around cylinder center (this must include initial conditions), the azimuth and axial components can be moved outside the time and radial integration and lead to a common ``2\pi{}z`` factor in both sides of the equation, which cancels out.

```math
\int_{0}^{z}\int_{0}^{2\pi}d\theta{}dz=2\pi{}z
```

The integration over radial coordinate introduces the ``rdr`` factor from the differential volume and we get the final form of the equation to integrate.

```math
\int_{s}^{n}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}rdtdr=
\int_{0}^{\tau}\int_{s}^{n}
\frac{\partial}{\partial{}r}
\left(rk\frac{\partial{}T}{\partial{}r}\right)drdt
```

Effecting the inner integration and moving out constant terms from the integrals we have

```math
\rho{}c_{p}\left(T_P^{\tau}-T_P^{0}\right)\int_{s}^{n}rdr=
\int_{0}^{\tau}
\left(rk\frac{\partial{}T}{\partial{}r}\right)\bigg\vert_{s}^{n}dt
```

Expanding the evaluation of the definite integral between control volume boundaries ``s`` and ``n`` and performing a Crank-Nicolson integration of the right-hand side one gets

```math
\begin{aligned}
    \frac{\rho{}c_{p}}{\tau}
    \left(T_P^{\tau}-T_P^{0}\right)
    \left(\frac{r_n^2}{2}-\frac{r_s^2}{2}\right)
    &=f\left[
    r_nk_n\frac{T_N^{\tau}-T_P^{\tau}}{\delta_{P,N}}-
    r_sk_s\frac{T_P^{\tau}-T_S^{\tau}}{\delta_{P,S}}
    \right]\\[8pt]
    &+(1-f)\left[
    r_nk_n\frac{T_N^{0}-T_P^{0}}{\delta_{P,N}}-
    r_sk_s\frac{T_P^{0}-T_S^{0}}{\delta_{P,S}}
    \right]
\end{aligned}
```

Some coefficients appearing in the above equations are now grouped. Notice that for thermal conductivity ``k`` which is a function of temperature, the corresponding time-step temperature must be used for its evaluation. For ``\beta_{j}`` the lower case ``j`` represents the evaluation at the interface with control volume ``J``, what is a very specific notation.

```math
\begin{aligned}
    \alpha_{P}  & = \frac{\rho{}c_{p}}{2\tau}\left(r_n^2-r_s^2\right)\\[8pt]
    \beta_{j}   & = \frac{r_jk_j}{\delta_{P,J}}
\end{aligned}
```

For conciseness we make ``g=(1-f)`` and simplify the expression with the new coefficients as

```math
-f\beta_{s}T_S+
(\alpha_{P}+f\beta_{n}+f\beta_{s})T_P
-f\beta_{n}T_N
=
g\beta_{s}T_S^{0}+
(\alpha_{P}-g\beta_{n}-g\beta_{s})T_P^{0}+
g\beta_{n}T_N^{0}
```

\subsection{Implicit implementation}

For the fully implicit time-stepping scheme ``f=1`` the expression reduces to

```math
-\beta_{s}T_S+
(\alpha_{P}+\beta_{n}+\beta_{s})T_P
-\beta_{n}T_N
=
\alpha_{P}T_P^{0}
```

where the following coefficients are identified

```math
\begin{aligned}
    a_{S} & = -\beta_{s}\\[8pt]
    a_{N} & = -\beta_{n}\\[8pt]
    a_{P} & = \alpha_{P}+\beta_{n}+\beta_{s}
\end{aligned}
```

and the standard format FVM discretization is reached

```math
a_ST_S + a_PT_P + a_NT_N = \alpha_{P}T_P^{0}
```

A condition for symmetry is that no flux traverses the center of the cylinder at ``r=0``. That implies that *south* derivatives in discrete form of the equation must vanish to enforce ``\dot{q}(0,t)=0``, so the first row of the problem is modified to

```math
a_1T_P + a_NT_N = \alpha_{P}T_P^{0}\quad\text{where}\quad{}a_1=\alpha_{P}+\beta_{n}
```

Over the external radius ``r=R`` a Robin boundary condition is imposed. In this case the heat flux ``\dot{q}=U(T_\infty-T_P)`` takes the place of *north* term in FVM discretization and the equation writes

```math
a_ST_S + a_RT_P = \alpha_{P}T_P^{0}+UT_\infty\quad\text{where}\quad{}a_R=\alpha_{P}+U+\beta_{s}
```

It must be noted here that ``U=Rh``, where the actual heat transfer coefficient is ``h``. This should be self-evident from a dimensional analysis.

### Spherical coordinates 1-D

In the case of spherical coordinates we start with a modification in divergence operator as follows

```math
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}=
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)
```

The expression is again integrated over time and the differential volume ``dV``.

```math
\int_{V}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}dtdV=
\int_{0}^{\tau}\int_{V}
\frac{1}{r^2}\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)dVdt
```

This differential volume is given by ``dV=r^2\sin\phi{}dr{}d\theta{}d\phi``. Since the problem is specified to be symmetric around sphere center (this must include initial conditions), the polar and azimuth components can be moved outside the time and radial integration and lead to a common ``4\pi`` factor in both sides of the equation, which cancels out.

```math
\int_{0}^{\pi}\int_{0}^{2\pi}\sin\phi{}d\theta{}d\phi=4\pi
```

The integration over radial coordinate introduces the ``r^2dr`` factor from the differential volume and we get the final form of the equation to integrate.

```math
\int_{s}^{n}\int_{0}^{\tau}
\rho{}c_{p}\frac{\partial{}T}{\partial{}t}r^2dtdr=
\int_{0}^{\tau}\int_{s}^{n}
\frac{\partial}{\partial{}r}
\left(r^2k\frac{\partial{}T}{\partial{}r}\right)drdt
```

After effecting the inner integration and moving out constant terms from the integrals and expanding the evaluation of the definite integral between control volume boundaries ``s`` and ``n`` and performing a Crank-Nicolson integration of the right-hand side one gets

```math
\begin{aligned}
    \frac{\rho{}c_{p}}{\tau}
    \left(T_P^{\tau}-T_P^{0}\right)
    \left(\frac{r_n^3}{3}-\frac{r_s^3}{3}\right)
    &=f\left[
    r_n^2k_n\frac{T_N^{\tau}-T_P^{\tau}}{\delta_{P,N}}-
    r_s^2k_s\frac{T_P^{\tau}-T_S^{\tau}}{\delta_{P,S}}
    \right]\\[8pt]
    &+(1-f)\left[
    r_n^2k_n\frac{T_N^{0}-T_P^{0}}{\delta_{P,N}}-
    r_s^2k_s\frac{T_P^{0}-T_S^{0}}{\delta_{P,S}}
    \right]
\end{aligned}
```

Some coefficients appearing in the above equations are now grouped. Notice that for thermal conductivity ``k`` which is a function of temperature, the corresponding time-step temperature must be used for its evaluation. For ``\beta_{j}`` the lower case ``j`` represents the evaluation at the interface with control volume ``J``, what is a very specific notation.

```math
\begin{aligned}
    \alpha_{P}  & = \frac{\rho{}c_{p}}{3\tau}\left(r_n^3-r_s^3\right)\\[8pt]
    \beta_{j}   & = \frac{r_j^2k_j}{\delta_{P,J}}
\end{aligned}
```

For conciseness we make ``g=(1-f)`` and simplify the expression with the new coefficients as

```math
-f\beta_{s}T_S+
(\alpha_{P}+f\beta_{n}+f\beta_{s})T_P
-f\beta_{n}T_N
=
g\beta_{s}T_S^{0}+
(\alpha_{P}-g\beta_{n}-g\beta_{s})T_P^{0}+
g\beta_{n}T_N^{0}
```

\subsection{Implicit implementation}

For the fully implicit time-stepping scheme ``f=1`` the expression reduces to

```math
-\beta_{s}T_S+
(\alpha_{P}+\beta_{n}+\beta_{s})T_P
-\beta_{n}T_N
=
\alpha_{P}T_P^{0}
```

where the following coefficients are identified

```math
\begin{aligned}
    a_{S} & = -\beta_{s}\\[8pt]
    a_{N} & = -\beta_{n}\\[8pt]
    a_{P} & = \alpha_{P}+\beta_{n}+\beta_{s}
\end{aligned}
```

and the standard format FVM discretization is reached

```math
a_ST_S + a_PT_P + a_NT_N = \alpha_{P}T_P^{0}
```

A condition for symmetry is that no flux traverses the center of the sphere at ``r=0``. That implies that *south* derivatives in discrete form of the equation must vanish to enforce ``\dot{q}(0,t)=0``, so the first row of the problem is modified to

```math
a_1T_P + a_NT_N = \alpha_{P}T_P^{0}\quad\text{where}\quad{}a_1=\alpha_{P}+\beta_{n}
```

Over the external radius ``r=R`` a Robin boundary condition is imposed. In this case the heat flux ``\dot{q}=U(T_\infty-T_P)`` takes the place of *north* term in FVM discretization and the equation writes

```math
a_ST_S + a_RT_P = \alpha_{P}T_P^{0}+UT_\infty\quad\text{where}\quad{}a_R=\alpha_{P}+U+\beta_{s}
```

It must be noted here that ``U=R^2h``, where the actual heat transfer coefficient is ``h``. This should be self-evident from a dimensional analysis.

## Mass transfer

### Carbon diffusion in plain iron

```math
\frac{\partial{}x}{\partial{}t}=\nabla\cdotp{}(D(x)\nabla{}T)
```

```math
\frac{\partial{}x}{\partial{}t}=
\frac{\partial}{\partial{}x}
\left(D(x)\frac{\partial{}x}{\partial{}z}\right)
```

```math
\int_{s}^{n}\int_{0}^{\tau}
\frac{\partial{}x}{\partial{}t}dtdz=
\int_{0}^{\tau}\int_{s}^{n}
\frac{\partial}{\partial{}z}
\left(D(x)\frac{\partial{}x}{\partial{}z}\right)dzdt
```

```math
\left(x_P^{\tau}-x_P^{0}\right)(w_{n}-w_{s})=
\int_{0}^{\tau}
\left(D(x)\frac{\partial{}T}{\partial{}r}\right)\bigg\vert_{s}^{n}dt
```

```math
\begin{align}
\left(x_P^{\tau}-x_P^{0}\right)\frac{(w_{n}-w_{s})}{\tau}&=
f\left[
D(x_n)\frac{x_N^{\tau}-x_P^{\tau}}{\delta_{P,N}}-
D(x_s)\frac{x_P^{\tau}-x_S^{\tau}}{\delta_{P,S}}
\right]\\[8pt]
&+(1-f)\left[
D(x_n)\frac{x_N^{0}-x_P^{0}}{\delta_{P,N}}-
D(x_s)\frac{x_P^{0}-x_S^{0}}{\delta_{P,S}}
\right]
\end{align}
```

```math
\begin{align}
\alpha_{P}  & = \frac{(w_{n}-w_{s})}{\tau}\\[8pt]
\beta_{j}   & = \frac{D(x_j)}{\delta_{P,J}}
\end{align}
```

```math
-f\beta_{s}x_S+
(\alpha_{P}+f\beta_{n}+f\beta_{s})x_P
-f\beta_{n}x_N
=
g\beta_{s}x_S^{0}+
(\alpha_{P}-g\beta_{n}-g\beta_{s})x_P^{0}+
g\beta_{n}x_N^{0}
```

```math
-\beta_{s}x_S+
(\alpha_{P}+\beta_{n}+\beta_{s})x_P
-\beta_{n}x_N
=
\alpha_{P}x_P^{0}
```

```math
\begin{align}
a_{S} & = -\beta_{s}\\[8pt]
a_{N} & = -\beta_{n}\\[8pt]
a_{P} & = \alpha_{P}+\beta_{n}+\beta_{s}
\end{align}
```

```math
a_Sx_S + a_Px_P + a_Nx_N = \alpha_{P}x_P^{0}
```

```math
a_1x_P + a_Nx_N = \alpha_{P}x_P^{0}\quad\text{where}\quad{}a_1=\alpha_{P}+\beta_{n}
```

```math
a_Sx_S + a_Rx_P = \alpha_{P}x_P^{0}+hx_\infty\quad\text{where}\quad{}a_R=\alpha_{P}+h+\beta_{s}
```

!!! note "About mass intake calculation"

```math
\rho_{Fe} = \frac{m_{Fe}}{V_{cell}}
```

```math
y_{C} = \frac{m_{C}}{m_{Fe} + m_{C}}
```

```math
m_{Fe+C} = \frac{m_{Fe}}{1 - y_{C}}
```

```math
\rho_{Fe+C} = \rho_{Fe}\frac{1}{1 - y_{C}}
```

```math
\sigma = \int_{0}^{L}\rho(z)y_{C}(z)dz
```

```math
\sigma = \rho_{Fe}\int_{0}^{L}\frac{y_{C}(z)}{1-y_{C}(z)}dz
```

```math
\Delta\sigma = \rho_{Fe}\left(\int_{0}^{L}\frac{y_{C}(z)}{1-y_{C}(z)}dz\right)\biggr\vert_{t=0}^{t=t_{f}}
```
