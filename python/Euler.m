function soln = Euler(h, x0, y0, interval_length, func)
nsteps = floor(interval_length/h) + 1;
x = zeros(nsteps,1);
y = zeros(nsteps,1);
x(1) = x0;
y(1) = y0;
for i=2:nsteps
y(i) = y(i-1) + h* func(x(i-1), y(i-1));
x(i) = x(i-1) + h;
end
soln = y(nsteps);