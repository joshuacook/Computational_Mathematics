function soln = runge_kutte(h, x0, y0, interval_length, func)
nsteps = floor(interval_length/h) + 1;
x = zeros(nsteps,1);
y = zeros(nsteps,1);
x(1) = x0;
y(1) = y0;
for i=2:nsteps
    k1(i) = func(x(i-1), y(i-1));
    k2(i) = func(x(i-1)+0.5*h, y(i-1)+0.5*h*k1(i));
    k3(i) = func(x(i-1)+0.5*h, y(i-1)+0.5*h*k2(i));
    k4(i) = func(x(i-1)+h, y(i-1)+h*k3(i));
y(i) = y(i-1) + h/6*(k1(i)+2*k2(i)+2*k3(i)+k4(i));
x(i) = x(i-1) + h;
end
soln = y(nsteps);