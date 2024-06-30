using MKL
using Statistics
using Plots
using LinearAlgebra
using Interpolations
using DataFrames
using GLM
using SparseArrays
#using Polynomials


lmda_values=[0.001, 0.01, 0.1, 1, 10]
c_star_vect=[0.22580876946449277, 0.37393441796302795, 0.775405079126358, 2.0734446942806244, 6.350513994693756]
#[0.22580876946449277, 0.37393441796302795, 0.775405079126358, 2.0734446942806244, 6.350513994693756]
#[0.2258, 0.3739, 0.7754, 2.073, 6.351]
N_vect=[1000, 1000, 2000, 3000, 4000]

K_vect=[]
K1_extracted_vect=[]
K_1_VECT=[]
K_2_VECT=[]
INTERCEPT_VECT=[]
F_eta_vect=[]
speed_vect=[]
#u_mean_vect=zeros((sum(N_vect)))


for i in 1:length(lmda_values)
    # Important Variables
    N=N_vect[i]
    one_vect=ones(N)
    lmda=lmda_values[i]
    C_star=c_star_vect[i]
    p=1
    h=0.0001
    Wave_Threshold=0.5
    u_vect=ones((N))
    Threshold_attained_time_series=zeros((N))
    # Define and initialise line lattice in abstract LDS framework
    global u_vect=zeros((N))
    initial_conditions=50
    for m in 1:initial_conditions
        u_vect[m]=1
    end

    # Defining matrix A (note, A will be of size NxN)
    d=-2*ones((N))
    v1=ones((N))
    v2=ones((N))
    A=Tridiagonal(v1[1:N-1], d, v2[1:N-1])
    A[1,2]=2
    A[N,N-1]=2
    global A=sparse(A)
    # Defining Runge-Kutta 4th order function
    function Runge_Kutta_4th_vect(mat_product::Array{Float64,1}, u_vect::Array{Float64,1})
        k_1=lmda*(mat_product)+(u_vect.^(p)).*(one_vect-u_vect)
        k_2=lmda*((mat_product+h*(k_1/2)))+((u_vect+h*(k_1/2)).^(p)).*(one_vect-(u_vect+h*(k_1/2)))
        k_3=lmda*((mat_product+h*(k_2/2)))+((u_vect+h*(k_2/2)).^(p)).*(one_vect-(u_vect+h*(k_2/2)))
        k_4=lmda*((mat_product+h*(k_3)))+((u_vect+h*(k_3)).^(p)).*(one_vect-(u_vect+h*(k_3)))
        u_vect=u_vect+(h/6)*(k_1+2*k_2+2*k_3+k_4)
        return u_vect
    end


    # Main while loop
    global t=1
    global wavesum=50
    while u_vect[N-100]<Wave_Threshold
        global t=t+1
        previous_wavesum=wavesum
        #print(" ")
        #print(u_vect[N-100])
        #print(" ")
        
        # Do matrix multiplication
        mat_product=A*u_vect
        # Now do rk_vect method
        global u_vect=Runge_Kutta_4th_vect(mat_product, u_vect)
        #global u_vect=sparse(u_vect)
        # Now record where wave passes threshold value in u_vect
        #print(" ")
        global wavesum=sum(u_vect .> 0.5)
        #print(wavesum)
        #print(" ")
        if (wavesum>previous_wavesum)
            Threshold_attained_time_series[wavesum]=h*t
            print(" ")
            print("Wave at node:"*string(wavesum))
        end

        

        # calculates u mean
        u_mean=mean(u_vect)

        
    end


    # Now calculating wavespeed
    Speed_series=zeros((N))
    for j in (initial_conditions+1):N
        # Boundary
        if (j==1 && Threshold_attained_time_series[j]!=0)
            Speed_series[j]=1/abs(Threshold_attained_time_series[j]-Threshold_attained_time_series[j+1])
        end
        if (j==N && Threshold_attained_time_series[j,1]!=0)
            Speed_series[j]=1/abs(Threshold_attained_time_series[j]-Threshold_attained_time_series[j-1])
        end
        # Interior
        if (j!=1 && j!=N && Threshold_attained_time_series[j]!=0)
            Speed_series[j]=1/abs(Threshold_attained_time_series[j+1]-Threshold_attained_time_series[j])
        end
    end

    append!(speed_vect, Speed_series[N-110])
    print(speed_vect)

    # C Vs. 1/t

    # Linear approximation 
    y=sort(Speed_series[500:N-110], rev=true)
    x=sort!(one_vect[500:N-110]./Threshold_attained_time_series[500:N-110])
    print(length(x))
    print(" stuff ")
    print(length(y))
    #DATAFRAME OF POINTS
    columns = Any[x, y];
    Speed_time_dataframe=DataFrame(columns, [:col1, :col2])
    print(Speed_time_dataframe)
    print(unique(y))
    uniq_df=unique(Speed_time_dataframe, :col2)
    print(uniq_df)
    filtered_time=uniq_df[:,1]
    filtered_speeds=uniq_df[:,2]
    itp=interpolate((filtered_time,), filtered_speeds, Gridded(Linear()))
  

    # Approximation of gradient K1
    data = DataFrame(X=filtered_time, Y=filtered_speeds)
    line = lm(@formula(Y ~ X), data)
    print(line)
    K_1= GLM.coef(line)[2] # EXTRACTING K1 HERE
    append!(K_1_VECT,K_1)
    intercept=GLM.coef(line)[1]
    append!(INTERCEPT_VECT, intercept) 

    print(" K1 IS: "*string(K_1))
    print(" INTERCEPT=c* IS: "*string(intercept))
    sle=GLM.coef(line)[2]*(one_vect[70:N-110]./Threshold_attained_time_series[70:N-110]).+GLM.coef(line)[1] 
    #GLM.coef(line)[2]*filtered_time.+GLM.coef(line)[1] 

    # C Vs. t 

    if i==1
        display(plot(Threshold_attained_time_series[1:N-110], Speed_series[1:N-110], label=["c"],xlabel="t",ylabel="c",legend=false))
        hline!([C_star], linestyle=:dash, label=["c*="*string(C_star)])
        #savefig(hline!([C_star], linestyle=:dash, label=["c*="*string(C_star)]), "combi_wkbj_pulled_speed_convergence"*string(lmda)*"_1000.png")
    end

    if i>1
        display(plot!(Threshold_attained_time_series[1:N-110], Speed_series[1:N-110], label=["c"],xlabel="t",ylabel="c",legend=false))
        #hline!([C_star], linestyle=:dash, label=["c*="*string(C_star)], legend=:topright)
        savefig(hline!([C_star], linestyle=:dash, label=["c*="*string(C_star)]), "combi_wkbj_pulled_speed_convergence"*string(lmda)*"_1000.png")
    end
    #display(plot(Threshold_attained_time_series[1:N-110], Speed_series[1:N-110],label=["c"],xlabel="t",ylabel="c",legend=:bottomright))
    #savefig(hline!([C_star], linestyle=:dash, label=["c*="*string(C_star)]), "wkbj_pulled_speed_convergence"*string(lmda)*"_1000.png")

    # C Vs. 1/t

    #display(plot(one_vect[1:N-110]./Threshold_attained_time_series[1:N-110], Speed_series[1:N-110],label=["c"],xlabel="1/t",ylabel="c",legend=:bottomright))
    #display(hline!([C_star], linestyle=:dash, label=["η="*string(C_star)]))
    

    #display(plot!(one_vect[70:N-110]./Threshold_attained_time_series[70:N-110],sle,linestyle=:dash, label=["c~K1*1/t+c*"]))
    #savefig(hline!([C_star], linestyle=:dash, label=["c*="*string(C_star)]), "wkbj_straight_line_speed_approximation"*string(lmda)*"_1000.png")

    # Approximation of K2
    #K_2 ~ t^3/2 * (c - η - K_1*1/t)
    K_2=(Threshold_attained_time_series[N-110])^(3/2)*(Speed_series[N-110] - C_star -K_1*(1/Threshold_attained_time_series[N-110]))
    append!(K_2_VECT, K_2)

    # Calculate F_eta:
    F_eta=asinh(Speed_series[N-110]/(2*lmda))
    append!(F_eta_vect, F_eta)

    

    
end

# Calculating F_eta*K_1
F_eta_K_1=K_1_VECT.*F_eta_vect
print(" END SPEEDs are: ")
print(speed_vect)
print(" K_1 VECT IS: ")
print(K_1_VECT)
print("   ")
print(" INTERCEPTS/η ESTIMATES ARE: ")
print(INTERCEPT_VECT)
print("   ")
print(" K_2 VECT IS: ")
print(K_2_VECT)
print("   ")
print(" APPROXIMATED F_ETA VECT IS: ")
print(F_eta_vect)
print("   ")
print(" APPROXIMATED Fη*K_1 VECT: ")
print(F_eta_K_1)
