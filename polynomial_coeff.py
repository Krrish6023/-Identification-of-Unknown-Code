import sympy
from sympy import symbols, expand, binomial

# ---------------------------------------------------------------------
# 1) Define Symbols and Parameters
# ---------------------------------------------------------------------
X, Y = symbols('X Y', complex=True)
# You can treat r, gamma, n, p, rho, N as either symbols or numeric values.
# If they are known numeric parameters, just set them as integers/floats.
r_val    = 2       # Example integer for r
gamma_val= 3       # Example integer for gamma
n_val    = 4       # Example integer for n
p_val    = 1       # Example integer for p
rho_val  = 2       # Example integer for rho
N_val    = 5       # Example integer for N

# If you want them as Sympy symbols, do:
# r, gamma, n, p, rho, N = symbols('r gamma n p rho N', positive=True)

# ---------------------------------------------------------------------
# 2) Define Q(X, Y) and P(Y)
#    From the image, a simplified reading is:
#    P(Y) = [ (1 + Y)^(rho^r) ]^n = (1 + Y)^(n * rho^r)
#    Q(X, Y) = Product_{i=1 to r} [ (1/2)*((1+Y)(1+X) + (1-Y)(1-X)) ]^(gamma^n)
#    But note that ((1+Y)(1+X) + (1-Y)(1-X))/2 = (1 + X*Y).
#    So effectively Q(X, Y) = (1 + X*Y)^(r * gamma^n).
# ---------------------------------------------------------------------
# We’ll just substitute in the numeric values. 
# If you prefer a purely symbolic approach for r, gamma, n, etc., see below.

Q_expr = (1 + X*Y)**(r_val * (gamma_val**n_val))
P_expr = (1 + Y)**(n_val * (rho_val**r_val))

# ---------------------------------------------------------------------
# 3) Function to get coefficients:
#    For a bivariate polynomial expr in X and Y, we want the coefficient 
#    of X^p * Y^t. Sympy can do this by repeated .coeff calls.
# ---------------------------------------------------------------------
def bivariate_coefficient(expr, x_var, x_pow, y_var, y_pow):
    """
    Extract the coefficient of x_var^x_pow * y_var^y_pow in expr.
    This uses Sympy's .coeff(...) repeatedly.
    """
    # Expand the polynomial (important for correct .coeff() extraction)
    expanded_expr = expand(expr)
    # First get the coefficient wrt x_var^x_pow, then wrt y_var^y_pow
    coeff_x_part = expanded_expr.coeff(x_var, x_pow)
    coeff_final  = coeff_x_part.coeff(y_var, y_pow)
    return coeff_final

# ---------------------------------------------------------------------
# 4) Compute the sum for a_w
#    a_w <= sum_{t=0 to N} [coef(Q, X^p Y^t) * coef(P, Y^t) / binomial(N, t)]
# ---------------------------------------------------------------------
a_w = 0
for t in range(N_val + 1):
    # coefficient of X^p_val Y^t in Q_expr
    coef_Q = bivariate_coefficient(Q_expr, X, p_val, Y, t)
    # coefficient of Y^t in P_expr
    # For P_expr, we only need a univariate coefficient extraction in Y.
    coef_P = expand(P_expr).coeff(Y, t)
    
    # Multiply by binomial(N_val, t) and accumulate
    a_w += coef_Q * coef_P / binomial(N_val, t)

print("Computed sum (upper bound for a_w):", a_w)
