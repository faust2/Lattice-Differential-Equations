using Statistics
using Plots

# Define and initialise line lattice
# Important Variables
global N=500
global lmda_values=[1]
global p=1
global h=0.01
global time=200
global runtime=floor(time/h)
global Wave_Threshold=0.5
global speed_averages=Vector{Float64}()
global multiple=1000

# Defining Runge-Kutta 4th order function
function Runge_Kutta_4th(cent,left,right,lmda,p)
    u=cent
    k_1=lmda*(right+left-2*u)+(u^(p))*(1-u)
    k_2=lmda*(right+left-2*(u+h*(k_1/2)))+((u+h*(k_1/2))^(p))*(1-(u+h*(k_1/2)))
    k_3=lmda*(right+left-2*(u+h*(k_2/2)))+((u+h*(k_2/2))^(p))*(1-(u+h*(k_2/2)))
    k_4=lmda*(right+left-2*(u+h*(k_3)))+((u+h*(k_3))^(p))*(1-(u+h*(k_3)))
    u=u+(h/6)*(k_1+2*k_2+2*k_3+k_4)
    return u
end

for r in 1:length(lmda_values)
    # Initialise Line Lattice
    
    global lmda=lmda_values[r]
    global Line_Lattice=zeros((N, 1))
    global Threshold_attained_time_series=zeros((N, 1))

    lower_bound=Int(N/2)-1
    upper_bound=lower_bound+1
    cent_coord=Int(N/2)

    global initial_conditions=50

    for m in 1:initial_conditions
        Line_Lattice[m,1]=1
    end
    print(" ")
    print("Lattice at coupling strength: "*string(lmda))

    # Main for loop
    for t in 1:runtime
        if (mod(t,multiple)==0)
            print(" ")
            print(t*h)
        end
        #print(t)
        #print(" ")
        global New_Lattice=zeros((N,1))

        for i in 1:N
            # Boundary Conditions
            # Left Boundary
            if (i==1)
                rk_output=Runge_Kutta_4th(Line_Lattice[i,1],Line_Lattice[i+1,1],Line_Lattice[i+1,1],lmda,p)
                New_Lattice[i,1]=rk_output
                if (New_Lattice[i,1]>=Wave_Threshold)
                    if (Threshold_attained_time_series[i,1]==0)
                        Threshold_attained_time_series[i,1]=h*t
                    end
                end  
            end

            # Right Boundary
            if (i==N)
                rk_output=Runge_Kutta_4th(Line_Lattice[i,1],Line_Lattice[i-1,1],Line_Lattice[i-1,1],lmda,p)
                New_Lattice[i,1]=rk_output
                if (New_Lattice[i,1]>=Wave_Threshold)
                    if (Threshold_attained_time_series[i,1]==0)
                        Threshold_attained_time_series[i,1]=h*t
                    end
                end 
            end

            # Domain interior
            if (i>1 && i<N)
                rk_output=Runge_Kutta_4th(Line_Lattice[i,1],Line_Lattice[i-1,1],Line_Lattice[i+1,1],lmda,p)
                New_Lattice[i,1]=rk_output
                if (New_Lattice[i,1]>=Wave_Threshold)
                    if (Threshold_attained_time_series[i,1]==0)
                        Threshold_attained_time_series[i,1]=h*t
                    end
                end
            end
        end

        global Line_Lattice=New_Lattice

        if mod(t, multiple)==0
            display(plot(1:N, Line_Lattice[:,1], linestyle=:dash, linewidth=5, title="μ="*string(lmda)*": Lattice at time: "*string(t*h)))
            savefig(plot(1:N, Line_Lattice[:,1], linestyle=:dash, linewidth=5, title="μ="*string(lmda)*": Lattice at time: "*string(t*h)), "Line Lattice at time_"*string(t*h)*".png")
        end
        
    end

    
    # Wavespeed computations
    # Speed=distance/time
    global Speed_series=zeros((N, 1))
    for j in (initial_conditions+1):N
        # Boundary
        if (j==1 && Threshold_attained_time_series[j,1]!=0)
            Speed_series[j,1]=1/abs(Threshold_attained_time_series[j,1]-Threshold_attained_time_series[j+1,1])
        end
        if (j==N && Threshold_attained_time_series[j,1]!=0)
            Speed_series[j,1]=1/abs(Threshold_attained_time_series[j,1]-Threshold_attained_time_series[j-1,1])
        end
        # Interior
        if (j!=1 && j!=N && Threshold_attained_time_series[j,1]!=0)
            Speed_series[j,1]=1/abs(Threshold_attained_time_series[j+1,1]-Threshold_attained_time_series[j,1])
        end
    end
    
    global average_speed=mean(Speed_series[:,1])
    push!(speed_averages, average_speed)

    #display(plot(1:N, Speed_series[:,1], legend=false, title="Lattice nodes wavespeeds at coupling strength: "*string(lmda)))
    #savefig(plot(1:4500, Speed_series[1:4500,1], lgened=false, xlabel="Lambda values", ylabel="Average Wavespeed"),"1d lattice lambda=1.png")
    display(plot(1:N, Threshold_attained_time_series[:,1], legend=false, title="Lattice nodes wavespeeds at coupling strength: "*string(lmda)))
    savefig(plot(1:4500, Threshold_attained_time_series[1:4500,1], lgened=false, xlabel="Lambda values", ylabel="Average Wavespeed"),"1d lattice lambda=1 time thresholds.png")
end
#display(plot(lmda_values, speed_averages, markershapes=:cross, xlabel="Lambda values", ylabel="Average Wavespeed"))
#savefig(plot(lmda_values, speed_averages, markershapes=:cross, xlabel="Lambda values", ylabel="Average Wavespeed"),"1d lattice wavespeed by lambda.png")
