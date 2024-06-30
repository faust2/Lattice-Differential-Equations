using MKL
using Statistics
using Plots
using LinearAlgebra
using CalculusWithJulia
import Contour: contours, levels, lines, coordinates
using ForwardDiff


mu_values=[0.001, 0.01, 0.1, 1, 10]
c_star_vect=[0.22580876946449277, 0.37393441796302795, 0.775405079126358, 2.0734446942806244, 6.350513994693756]


# Graph 1 c* vs mu
savefig(plot(mu_values, c_star_vect, xlabel="μ", ylabel="c*", legend=false, title="c* VS. μ"), "1pulled case, c_star by mu.png")


# Graph 2 lambda vs mu
#lambda_vec=[asinh(η/(2*μ))]
savefig(plot(mu_values, asinh.(c_star_vect./(2*mu_values)), xlabel="μ", ylabel="λ", legend=false, title="λ VS. μ"), "222pulled case, lambda by mu.png")


#f(η,μ)=asinh(η/(2*μ))
#surface(mu_values, eta_values, f, xlabel="μ", ylabel="η",zlabel="λ", camera = (70, 30))
#savefig(plot(surface(mu_values, eta_values, f, xlabel="μ", ylabel="η",zlabel="λ", camera = (10, 30))), "1d_sing_sol_2.png")
