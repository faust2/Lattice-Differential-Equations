using MDBM
using LinearAlgebra
using Plots


mu_vec=[1, 0.01, 0.0001, 0.000001]

for i in 1:length(mu_vec)
    # Define coupling strength
    global μ=mu_vec[i]
    function f(ξ, η)
        return μ*(exp(-asinh(ξ/(2*μ)))+exp(asinh(ξ/(2*μ)))+exp(-asinh(η/(2*μ)))+exp(asinh(η/(2*μ)))-4)+1-ξ*asinh(ξ/(2*μ))-η*asinh(η/(2*μ))
    end

    global ax1=Axis(collect(range(-2.1,2.1,1000)),"x") # initial grid in x direction
    global ax2=Axis(collect(range(-2.1,2.1,1000)),"y") # initial grid in y direction
    global mymdbm=MDBM_Problem(f,[ax1,ax2])
    global iteration=1 #number of refinements (resolution doubling)
    solve!(mymdbm, iteration)

    #interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
    global x_sol,y_sol=getinterpolatedsolution(mymdbm)

    # normalizing betwen -1 and 1
    #if i==1
        #global x_sol_normed_1 = (x_sol .- minimum(x_sol)) / (maximum(x_sol) - minimum(x_sol)).*(1+1).-1
        #global y_sol_normed_1 = (y_sol .- minimum(y_sol)) / (maximum(y_sol) - minimum(y_sol)).*(1+1).-1
    #end
    #if i>1
        #global total_x_sol=[]
        #global total_y_sol=[]
        #global x_sol_normed = (x_sol .- minimum(x_sol)) / (maximum(x_sol) - minimum(x_sol)).*(1+1).-1
        #global y_sol_normed = (y_sol .- minimum(y_sol)) / (maximum(y_sol) - minimum(y_sol)).*(1+1).-1
    #end

    global s=3
    global x_sol_normed = (x_sol .- minimum(x_sol)) / (maximum(x_sol) - minimum(x_sol)).*(1+1).-1
    global y_sol_normed = (y_sol .- minimum(y_sol)) / (maximum(y_sol) - minimum(y_sol)).*(1+1).-1
    if i==1
        
        display(scatter(x_sol_normed, y_sol_normed, xlabel="ξ", ylabel="η", markersize=s, legend = :outertopright, label="μ="*string(μ)))
    end
    if i>1
        #display(scatter!(x_sol_normed, y_sol_normed, xlabel="ξ", ylabel="η", markersize=s, legend = :outertopright, label="μ="*string(μ)))
        savefig(scatter!(x_sol_normed, y_sol_normed, xlabel="ξ", ylabel="η", markersize=s, legend = :outertopright, label="μ="*string(μ)), "2D Square numerics.png")
    end



end









