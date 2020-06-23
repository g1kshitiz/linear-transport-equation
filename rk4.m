% Implementation of 4th order Runge Kutta
function y = rk4(x, t, phi, dx, dt, c)
   y = phi;
   for n = 1:(length(t)-1)
       for i = 3:(length(x)-2)
            k1 = (y(i+1,n)- y(i-1,n))/(2*dx);
            k2 = (y(i+1,n) + y(i-1,n) - 2*y(i,n))/dx^2;
            k3 = (y(i+2,n) - y(i-2,n) - 2*y(i+1,n) + 2*y(i-1,n) )/(2*dx^3);
            k4 = (y(i+2,n) + y(i-2,n) - 4*y(i+1,n) -4*y(i-1,n) +y(i,n))/dx^4;
            % step wise solution for the pseudo-ode
            Rn = -c*k1;
            R1 = Rn + c^3*dt/3*k2;
            R2 = Rn + c^2*dt/2*k2 + c^3*dt^2/4*k3;
            R3 = Rn + c^3*dt*k2 - c^3 * dt^2*k3/2 - c^4 * dt^3 * k4/4;
            y(i, n+1) = y(i,n) + dt*(Rn + 2*R1 + 2*R2 + R3)/6;
       end
   end