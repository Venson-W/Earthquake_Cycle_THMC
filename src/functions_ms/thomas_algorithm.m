
% Thomas algorithm for solving tridiagonal system
% a: lower diagonal, b: main diagonal, c: upper diagonal, d: right-hand side
function x = thomas_algorithm(a, b, c, d, n)
    % Check inputs
    if n <= 0
        error('System size must be positive');
    end
    
    % Initialize solution vector
    x = zeros(n, 1);
    
    % Special case for n=1
    if n == 1
        x(1) = d(1) / b(1);
        return;
    end
    
    % Forward elimination
    c_prime = zeros(n, 1);
    d_prime = zeros(n, 1);
    
    c_prime(1) = c(1) / b(1);
    d_prime(1) = d(1) / b(1);
    
    for i = 2:n
        denominator = b(i) - a(i) * c_prime(i-1);
        if abs(denominator) < 1e-10
            denominator = 1e-10 * sign(denominator);  % Avoid division by zero
        end
        c_prime(i) = c(i) / denominator;
        d_prime(i) = (d(i) - a(i) * d_prime(i-1)) / denominator;
    end
    
    % Back substitution
    x(n) = d_prime(n);
    
    for i = n-1:-1:1
        x(i) = d_prime(i) - c_prime(i) * x(i+1);
    end
