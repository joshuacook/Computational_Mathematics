<span>**MATH 481A Numerical Analysis   Fall 2014**</span>

<span>**Homework and Lab Assignment \# 5 – ODEs Due: Tue., Dec. 16,
2014**</span>

In this assignment you will write a Python script to solve ODEs using a
second order Runge-Kutta method. You will also analyze and verify the
order of the method, and use your code to solve a first order ODE and
the linearized pendulum equation, a second order ODE, by writing it as a
first order system.

You can work with others and discuss the problems, but each student must
write his/her own, independent solution. If you are unsure about what i
mean by this, please ask!

<span>**Second order Runge–Kutta Method:**</span> Given the first order
ODE [eq:ODE] = f(y(t),t),    t (t~0~, t~f~],    y(t~0~) = y~0~ we
consider the two-stage Runge-Kutta method:

[eq:RK2]

$$\begin{aligned}
k_1 \ =  & \, \Delta t f(y_n,t_n) \label{eq:RKk1} \\[0.1in] 
k_2 \ = & \, \Delta t f(y_n + \frac{k_1}{2},t_n + \frac{\Delta t}{2}) \label{eq:RKk2} \\[0.1in] 
y_{n+1} = & \ y_n + k_2. \label{eq:RKmp}\end{aligned}$$

to approximate its solution.

<span>**Problem 1.**</span> Find the local truncation error of the
Runge-Kutta method by comparing the Taylor expansion of $y(t_{n+1})$ and
that of
$y(t_n) + \Delta t f(y(t_n) + \frac{k_1}{2}, t_n + \frac{\Delta t}{2})$,
around $t_n$, noting that the Taylor expansion of $f(y + h, t + k)$
around $(y,t)$ is

$$\begin{aligned}
f(y + h, t + k) = & f(y, t) + h \frac{\partial f}{\partial y}(y,t) + k \frac{\partial f}{\partial t}(y,t) \no \\[0.1in]
 + & \half \left( h^2 \frac{\partial^2 f}{\partial y^2}(y,t) + 2 h k \frac{\partial^2 f}{\partial y \partial t}(y,t) + k^2 \frac{\partial^2 f}{\partial t^2}(y,t) \right) + {\cal O}(h^3,k^3). \no\end{aligned}$$

<span>**Problem 2.**</span> Write a function, `RK2_step`, that
implements one step of the above Runge-Kutta method. The function should
take as input the name of the function $f(y(t), t)$ describing the right
hand side of equation , $y_n$, $t_n$, and $\Delta t$, and returns
$y_{n+1}$. <span>**Problem 3.**</span> Write a function called
`RK2_method` within the same file as the function in problem \#2) that
takes as input the name of the python function defining $f(y(t),t)$,
$t_0$, $t_f$, $y_0$, and $\Delta t$ and returns as output two 1d `numpy`
arrays `y` and `t` containing the approximate solution of the ODE and
the time values it was computed at. Your function should do the
following:

Determine the number, $N$, of values $t_n \in [t_0, t_f]$ at which the
solution will be calculated, calculate those values, and store them in
the 1d array, `t`. These values are $t_n = t_0 + n \Delta t$ for
$n = 0, 1, 2, \dots, N-1$ and $T_N = t_f$. You should make sure the last
time step size passed to the function `RK2_step`, $\Delta t^*$, is
calculated so that $t_f = t_{N-1} + \Delta t^*$.

Create the 1d array of size $N$, `y`, to hold the solution at the
corresponding values stored in `t`.

Call the function `RK2_step` $N$ times to calculate $y_{n+1}$ for
$n = 0, 1, \dots N$. <span>**Problem 4.**</span> In the same file as you
wrote the functions in problems \#1 and \#2, write a function describing
$\displaystyle f(y(t),t) = 4 \, y\, (1 - y)$ and add the necessary code
to solve with this $f(y,t)$ as the right hands side, $t_0 = 0$,
$t_f = 1.0$, and $y_0 = 0.1$ using those functions, and plot the
computed solutions as needed to obtain:

Plots of the solution ($y \, vs. \, t$) calculated with
$\Delta t = 0.125, 0.0625, 0.03125$.

An estimate of the order of the method using Aitken’s extrapolation.

<span>**Problem 5.**</span> Consider the linearized pendulum equation
[eq:] = - ,    (0) = ~0~,   ’(0) = v~0~ where $g = 9.81 m/s^2$ stands
for the acceleration of gravity, $L$ is the length of the rod holding
the pendulum’s bob, and theta is the angle displacement measured
counter-clockwise from the negative $y$ axis.

$(i)$ Verify that $\theta(t) = a \cos{\omega t} + b \sin{\omega t}$,
where $\omega = g/L$ satisfies the equation, and $(ii)$ determine $a$
and $b$ in terms of the initial conditions, and $(iii)$ find the period,
$T$, of the oscillations of the pendulum (<span>*i.e.*</span>, find $T$
such that $\theta(t + T) = \theta(t)$).

Write the 2nd order ODE as a first order system of ODEs.

Write a function, `pendulum_exact` that takes as input $L$, $t_0$, $T$,
$\theta_0$, $v_0$, and $\Delta t$, and returns the exact solution of the
pendulum equation computed at the times $t_n = t_0 + n \Delta t$ for
$t \in [0, 2T]$.

Write a python function named `fpendulum` describing the right hand side
of your system and solve the equation using `RK2_method` for
$t \in [0, 2T]$ for the following conditions (in all cases use
$\Delta t = 0.1$):

$L = 2$, $\theta(0) = \frac{\pi}{6}$, $v(0) = 0$.

$L = 1$, $\theta(0) = \frac{\pi}{3}$, $v(0) = 0$.

$L = 1$, $\theta(0) = 0$, $v(0) = \frac{\pi}{10}$

Compute the exact solution for the same conditions and for each case
produce two plots: $(i)$ $y_{exact} \, vs. \, t$ and
$y_{RK2} \, vs. \, t$ (in the same axis), and $(ii)$ the error
$|y_{RK2} - y_{exact}|$. note that in this case, `RK2_method` should
return the 1d array `t` and a 2d array of size $(N+1) \times 2$ holding
the values of $\theta'$ and $\theta$ in columns 1 and 2 respectively.
<span>**Submission:**</span> Upload a `.pdf` file with the solution of
these problems and a files containing the functions `RK2_step`,
`RK2_method`, `fpendulum`, and the additional code you wrote to solve
the pendulum equation to
[](http://moodle.csun.edu/mod/assign/view.php?id=1713170)

