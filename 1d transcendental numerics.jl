using Roots
using ForwardDiff
using Polynomials
using IntervalArithmetic, IntervalRootFinding
mu=0.001
f(x) = x*(log((x+sqrt(x^2+4*mu^2))/(2*mu)))-sqrt(x^2+4*mu^2)+2*mu-1
function bisection(f, a, b, tol)
    if sign(f(a))==sign(f(b))
        error("Incorrect boundary values, f(x) must be between positive and negative.")
    end
    midpoint=(a+b)/2
    while abs(f(midpoint))>tol
        sign(f(midpoint))==sign(f(a)) ? a=midpoint : b=midpoint
        midpoint=(a+b)/2
    end
    return midpoint
end


# mu=0.001
bisection(f, 0.2, 0.25, 0.00000001)
# mu=0.01
#bisection(f, 0.3, 0.4, 0.00000001)
# mu=0.1
#bisection(f, 0.5, 1, 0.00000001)
# mu=1
#bisection(f, 2, 3, 0.00000001)
# mu=10
#bisection(f, 6, 7, 0.00000001)