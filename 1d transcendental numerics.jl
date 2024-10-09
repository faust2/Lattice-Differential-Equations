using Roots
using Plots
using ForwardDiff
using Polynomials
using IntervalArithmetic, IntervalRootFinding
mu=10
#μ=1

#global η=0
#f(ξ)=-ξ*log(ξ/2+sqrt((ξ/(2*μ)^2)+1))+μ*(log(ξ/2+sqrt((ξ/(2*μ)^2)+1)))+μ*(log(ξ/2+sqrt((ξ/(2*μ)^2)+1)))^(-1)-4*μ+1+μ+μ
#f(ξ, η)=-ξ*log(ξ/2+sqrt((ξ/(2*μ)^2)+1))+μ*(log(ξ/2+sqrt((ξ/(2*μ)^2)+1)))+μ*(log(ξ/2+sqrt((ξ/(2*μ)^2)+1)))^(-1)+μ*(log(η/2+sqrt((η/(2*μ)^2)+1)))-η*(log(η/2+sqrt((η/(2*μ)^2)+1)))+μ*(log(η/2+sqrt((η/(2*μ)^2)+1)))^(-1)-4*μ+1


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

#print(f(10))
#print(" ")
#print(f(0.1))
#print(" ")
#bisection(f, 10, 0.1, 0.0000000001)

# mu=0.001
#bisection(f, 0.2, 0.25, 0.00000001)
# mu=0.01
#bisection(f, 0.3, 0.4, 0.00000001)
# mu=0.1
#bisection(f, 0.5, 1, 0.00000001)
# mu=1
#bisection(f, 2, 3, 0.00000001)
# mu=10
#bisection(f, 6, 7, 0.00000001)

# Loop to produce all the c* values
mus=collect(range(0.001,10,1000))
global sped_vec=[]
for i in 1:length(mus)
    print(" ")
    print(i)
    mu=mus[i]
    f(x) = x*(log((x+sqrt(x^2+4*mu^2))/(2*mu)))-sqrt(x^2+4*mu^2)+2*mu-1
    sped=bisection(f, 10, 0.1, 0.0000000001)
    global sped_vec=push!(sped_vec, sped)

end
# Graph 1 c* vs mu
display(plot(mus, sped_vec,  linewidth=5, thickness_scaling = 1, c=:red, xlabel="Γ", ylabel="κ", legend=false, title="κ VS. Γ"))
savefig(plot(mus, sped_vec,  linewidth=5, thickness_scaling = 1, c=:red, xlabel="Γ", ylabel="κ", legend=false, title="κ VS. Γ"), "pulled case, c_star by mu.png")
# Graph 2 lambda vs mu
#lambda_vec=[asinh(η/(2*μ))]
display(plot(mus, asinh.(sped_vec./(2*mus)),  linewidth=5, thickness_scaling = 1, c=:red, xlabel="Γ", ylabel="ρ", legend=false, title="ρ VS. Γ"))
savefig(plot(mus, asinh.(sped_vec./(2*mus)),  linewidth=5, thickness_scaling = 1, c=:red, xlabel="Γ", ylabel="ρ", legend=false, title="ρ VS. Γ"), "pulled case, lambda by mu.png")
