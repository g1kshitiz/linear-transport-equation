% Implementation of Lax-Wendroff Method
function y = lax_wendroff(x, t, phi, c_bar)
    y = phi;
    for n = 1:(length(t)-1)
        for i = 2:(length(x)-1)
            y(i, n+1) = y(i,n) - c_bar/2 * ( y(i+1,n)-y(i-1,n) ) + c_bar^2/2 * (y(i+1,n) - 2*y(i, n) + y(i-1, n) );
        end 
    end
end