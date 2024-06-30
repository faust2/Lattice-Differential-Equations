using MKL
using Statistics
using Plots
using LinearAlgebra

lmda_values=[1]
#[0.1, 1, 10]
#[0.001, 0.01, 0.1, 1, 10]
#c_star_vect=[6.350513994693756]
#[0.775405079126358, 2.0734446942806244, 6.350513994693756]
#n_prime_vect=[1000]
N_vect=[1000]
alpha_values=[0.1]

#collect(range(0.001,10,20))
#point_point_point_one=[0.001, 0.001, 0.001, 0.001, 0.001]
#point_point_one=[0.01,0.01,0.01,0.01,0.01]
#point_one=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
#onees=[1,1,1,1,1,1,1,1,1,1]
#tens=[10,10,10,10,10,10,10,10,10,10]
#[0.001, 0.01, 0.1, 1, 10]
#alpha_values=collect(range(0.0001,10,20))
#lmda_values=collect(range(0.0001,10,20))
K_vect=[]
F_eta_vect=[]
c_dagger_vect=[]


#u_mean_vect=zeros((sum(N_vect)))
for k in 1:length(lmda_values)
    lmda=lmda_values[k]
    speed_vect=[]
    #c_star=c_star_vect[k]
    
    for i in 1:length(alpha_values)
        # Important Variables
        N=N_vect[i]
        one_vect=ones(N)
        alpha=alpha_values[i]
        print(" PARAMETER CYCLE: "*string(k)*" alpha is: "*string(alpha)*" coupling strength is: "*string(lmda)*". Beginning new count now:  ")
        #C_star=c_star_vect[i]
        p=1
        h=0.00001
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
        # Defining Runge-Kutta 4th order function
        function Runge_Kutta_4th_vect(mat_product::Array{Float64,1}, u_vect::Array{Float64,1})
            k_1=lmda*(mat_product)+(u_vect.^(p)).*(one_vect-u_vect).*(one_vect+u_vect./alpha)
            k_2=lmda*((mat_product+h*(k_1/2)))+((u_vect+h*(k_1/2)).^(p)).*(one_vect-(u_vect+h*(k_1/2))).*(one_vect+(u_vect+h*(k_1/2))./alpha)
            k_3=lmda*((mat_product+h*(k_2/2)))+((u_vect+h*(k_2/2)).^(p)).*(one_vect-(u_vect+h*(k_2/2))).*(one_vect+(u_vect+h*(k_2/2))./alpha)
            k_4=lmda*((mat_product+h*(k_3)))+((u_vect+h*(k_3)).^(p)).*(one_vect-(u_vect+h*(k_3))).*(one_vect+(u_vect+h*(k_3))./alpha)
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
            #u_mean=mean(u_vect)

            
        end

        # Item 2

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
        append!(c_dagger_vect, (2*sqrt(lmda/alpha)/pi))
        print(" CALC SPEED IS: "*string(speed_vect))
        moving_index=length(c_dagger_vect)
        c_dagger=c_dagger_vect[moving_index]

        # Extract K:
        #print(" ")
        #print(" TIME OF PRELIM SPEED IS: ")
        #print(string(Threshold_attained_time_series[N-111]))
        #print(" ")
        #print(string(Speed_series[N-111]))
        #print(" ")
        #print(" PRELIM K IS: ")
        #print(string(-(1/(Threshold_attained_time_series[N-111]))*log(Speed_series[N-111]-0.5)))
        #print(" TIME OF FINAL SPEED IS: ")
        #print(string(Threshold_attained_time_series[N-110]))
        #print(" ")

        
        #print(log(Speed_series[N-110]-((sqrt((4*lmda/alpha)-1)/2)*(1/atan(sqrt((4*lmda/alpha)-1))))))
        #print(" ")
        #print(Speed_series[N-110]-((sqrt((4*lmda/alpha)-1)/2)*(1/atan(sqrt((4*lmda/alpha)-1)))))
        #print(" ")
        print(" ")
        print(" c is: ")
        print(string(Speed_series[N-105]))
        print(" ")
        print(" c† is: ")
        print(string(Speed_series[N-110]))
        K=-(1/(Threshold_attained_time_series[N-110]))*log((Speed_series[N-110]-Speed_series[N-100]))
        print(" ")
        print(" K is: ")
        print(string(K))
        #(log(Speed_series[N-110]-(1/log(1/lmda)+(alpha*log(1/0.1)-log(alpha))/(log(1/0.1)^2))))
        #K=(1/(Threshold_attained_time_series[N-110]))*(log(c_star-Speed_series[N-110]))
        append!(K_vect, K)


        # Now plotting
        #display(plot(Threshold_attained_time_series[1:N-110], Speed_series[1:N-110],label=["C"],xlabel="Time",ylabel="C"))
        display(plot(Threshold_attained_time_series[100:N-105], (-one_vect[100:N-105]./Threshold_attained_time_series[100:N-105]).*log.((Speed_series[N-110].-(Speed_series[100:N-105]))),xlabel="t",ylabel="-(1/t)*ln(c†-c)",legend=:right))
        #display(hline!([(K)], linestyle=:dash, label=["K="*string(K)]))
        #c†
        #savefig(hline!(hline!([(K)], linestyle=:dash, label=["K="*string(K)])),"pushed_1OVERtln(c†-c)vt_"*string(lmda)*"alpha_value="*string(alpha)*"_.png")
        #s_vect=abs.(Speed_series.-c_dagger)
        #display(plot(Threshold_attained_time_series[1:N-110], s_vect[1:N-110], xlabel="Time",ylabel="(C†-C)",legend=false))
        #savefig(plot(Threshold_attained_time_series[1:N-110], s_vect[1:N-110], xlabel="Time",ylabel="(C†-C)",legend=false), "pushed_c_dag_minus_c"*string(lmda)*"alpha_value="*string(alpha)*".png")
        #mul_vect=abs.(one_vect./log.(Threshold_attained_time_series).*log.(s_vect))
        #display(plot(Threshold_attained_time_series[1:N-110], mul_vect[1:N-110], xlabel="Time",ylabel="1/log(t)(C†-C)", legend=false))
        #path="pushed_1log(t)(c_dag-c)_coupling_strength="*string(lmda)*"alpha_value="*string(alpha)*".png"
        #savefig(plot(Threshold_attained_time_series[1:N-110], mul_vect[1:N-110], xlabel="Time",ylabel="1/log(t)(C†-C)", legend=false),path)

       

        # Calculate F_eta:
        F_eta=asinh(c_dagger/(2*lmda))
        append!(F_eta_vect, F_eta)
    end

    #KF_eta=K_vect.*F_eta_vect
    print(" K_vect is: "*string(K_vect))
    #print(" Approximated F_eta is: "*string(F_eta_vect))
    #print(" K*Approximated F_eta is: "*string(KF_eta))

    #if (k==1)
        #display(plot(alpha_values, speed_vect, markershapes=:cross, xaxis=:log, yaxis=:log, xlabel="Alpha values", ylabel="C", legend=false))
    #end
    #ylim = [0, 0.3], yticks=0:0.02:0.3,
    #display(plot!(alpha_values, speed_vect, markershapes=:cross, xaxis=:log, yaxis=:log,  xlabel="Alpha values", ylabel="C", legend=false))
    #savefig(plot!(alpha_values, speed_vect, markershapes=:cross, xaxis=:log, yaxis=:log,  xlabel="Alpha values", ylabel="C", legend=false),"PushedPulled 1D Lattice wavespeed figure4b.png")

end


print(" K_vect is: "*string(K_vect))



# Plotting speed vs alpha, with limits set in 
#display(plot(alpha_values, c_dagger_vect[1:10], xlabel="alpha values", ylabel="C", xaxis=:log,  yaxis=:log, legend=false))
#plot!(alpha_values[1:6], 2*sqrt.(point_one[1:6]./alpha_values[1:6])/pi, linestyle=:dash) # this is the limit of alpha -> 0
#display(hline!([c_star_vect[1]], linestyle=:dash))

#plot!(alpha_values, c_dagger_vect[11:20], xlabel="alpha values", ylabel="C", xaxis=:log,  yaxis=:log, legend=false)
#plot!(alpha_values[1:6], 2*sqrt.(onees[1:6]./alpha_values[1:6])/pi, linestyle=:dash) # this is the limit of alpha -> 0
#display(hline!([c_star_vect[2]], linestyle=:dash))

#plot!(alpha_values, c_dagger_vect[21:30], xlabel="alpha values", ylabel="C", xaxis=:log,  yaxis=:log, legend=false)
#plot!(alpha_values[1:6], 2*sqrt.(tens[1:6]./alpha_values[1:6])/pi, linestyle=:dash) # this is the limit of alpha -> 0
#display(hline!([c_star_vect[3]], linestyle=:dash))

#savefig(hline!([c_star_vect[3]], linestyle=:dash),"alpha_wavespeed_mu.png")

#plot!(alpha_values, c_dagger_vect[16:20], xlabel="alpha values", ylabel="C", xaxis=:log,  yaxis=:log, legend=false)
#plot!(alpha_values, 2*sqrt.(onees./alpha_values)/pi, linestyle=:dash) # this is the limit of alpha -> 0
#display(hline!([c_star_vect[4]], linestyle=:dash))

#plot!(alpha_values, c_dagger_vect[21:25], xlabel="alpha values", ylabel="C", xaxis=:log,  yaxis=:log, legend=false)
#plot!(alpha_values, 2*sqrt.(tens./alpha_values)/pi, linestyle=:dash) # this is the limit of alpha -> 0
#display(hline!([c_star_vect[5]], linestyle=:dash))



#savefig(hline!([c_star_vect[5]], linestyle=:dash),"alpha_wavespeed_mu.png")
#display(hline!(2*sqrt.(c_star_vect./alpha_values)/pi))
# Plot K*F_eta Vs. Lmda:
KF_eta=K_vect.*F_eta_vect
print(" K_vect is: "*string(K_vect))
print(" Approximated F_eta is: "*string(F_eta_vect))
print(" K*Approximated F_eta is: "*string(KF_eta))
#display(plot(lmda_values, KF_eta, xlabel="Coupling Strength",ylabel="K(F_eta)",legend=false))
#savefig(plot(lmda_values, KF_eta, xlabel="Coupling Strength",ylabel="K(F_eta)",legend=false), "item 2 vers wiggles.png")