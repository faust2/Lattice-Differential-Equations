using MKL
using Statistics
using Plots
using LinearAlgebra
using SparseArrays

# Important Variables

N=200
percentage_stop=0.55
#lmda=1
lmda_vec=[0.001]
#c_star_vect=[0.22580876946449277, 0.37393441796302795, 0.775405079126358, 2.0734446942806244, 6.350513994693756]
p=1
h=0.001
Wave_Threshold=0.5
sorted_speed_by_angle_leading_point=[]

for z in 1:length(lmda_vec)
    global lmda=lmda_vec[z]
    print(" CYCLE:"*string(z)*" coupling strength="*string(lmda)*" ")
    p_θ_vect=[]
    global Threshold_attained_time_series=zeros(N^2)
    one_vect=ones((N^2))
    # Define and initialise lattice encoding vector, define stopping element
    init_lattice=spzeros(N,N)
    cent_coordinate=ceil(Int, N/2)
    cent_coord_x=1
    cent_coord_y=1
    #init_radius=ceil(Int, N*0.001)
    init_lattice[cent_coord_x, cent_coord_y]=1

    #for i in 1:50
        #init_lattice[i,1]=1
        #init_lattice[1,i]=1
    #end

    #for i in cent_coordinate-init_radius:cent_coordinate+init_radius
        #for j in cent_coordinate-init_radius:cent_coordinate+init_radius
            #init_lattice[i,j]=1
        #end
    #end

    print("check 1")

    # Heatmap of initial lattice
    he_1=heatmap(1:size(init_lattice,1), 1:size(init_lattice,2), init_lattice, aspect_ratio=:equal, xlims=(1, N), ylims=(1,N), c=cgrad(:viridis), title="Heatmap of Square Lattice Domain: Γ="*string(lmda))
    display(he_1)
    gui(he_1)
    savefig(he_1, "2DHeatmap_default_lambda"*string(lmda)*".png")

    print(" check 2")

    global u_vect=SparseVector{Float64, Int64}(vec(transpose(init_lattice)))

    # Stopping element

    global stopping_index=floor(Int, 30*N-(N/2))
    #floor(Int, 0.5*(N)^2-N)


    # Defining (sparse!) matrix A (note, A will be of size N^2xN^2)
    # To insert custom diagonal vectors, vec, use A[diagind(A, n)]=vec, where n=0 is central diagonal, 
    # and -n integers denote the diagonals to the left of the central diagonal, +n integers the right.

    # Matrix Shell
    A=spzeros(N^2,N^2)
    # Center diagonal
    v_c=ones(N^2)*-4
    A[diagind(A,0)]=v_c
    # Far right diagonal
    v_FR=ones(N^2-N)
    v_FR[1:N].=2
    # Right diagonal
    c_1=ones(N)
    c_1[1]=2
    c_1[N]=0
    v_R=repeat(c_1, outer = N)
    pop!(v_R) # gets rid of last element
    # Far left diagonal
    v_FL=reverse(v_FR)
    # Left diagonal
    v_L=reverse(v_R)
    # Creating true matrix A
    A[diagind(A,N)]=v_FR
    A[diagind(A,-N)]=v_FL
    A[diagind(A,1)]=v_R
    A[diagind(A,-1)]=v_L

    print(" check 3")

    #global mat_product=A*u_vect
    # Defining Runge-Kutta 4th order function
    # mat_product::Vector{Float64}, u_vect::Vector{Float64}
    # mat_product::SparseVector{Float64, Int64}, u_vect::Base.ReshapedArray{Float64, 1, Transpose{Float64, SparseMatrixCSC{Float64, Int64}}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}
    #u_vect::SparseVector{Float64, Int64}
    function Runge_Kutta_4th_vect(mat_product::SparseVector{Float64, Int64}, u_vect::SparseVector{Float64, Int64})
        k_1=lmda*(mat_product)+(u_vect.^(p)).*(one_vect-u_vect)
        k_2=lmda*((mat_product+h*(k_1/2)))+((u_vect+h*(k_1/2)).^(p)).*(one_vect-(u_vect+h*(k_1/2)))
        k_3=lmda*((mat_product+h*(k_2/2)))+((u_vect+h*(k_2/2)).^(p)).*(one_vect-(u_vect+h*(k_2/2)))
        k_4=lmda*((mat_product+h*(k_3)))+((u_vect+h*(k_3)).^(p)).*(one_vect-(u_vect+h*(k_3)))
        u_vect=u_vect+(h/6)*(k_1+2*k_2+2*k_3+k_4)
        return u_vect
    end

    #u_vect::Vector{Float64}
    #function fwd_eul_vect(mat_product::SparseVector{Float64, Int64}, u_vect::Vector{Float64})
        #u_vect=u_vect+h*(mat_product+(u_vect.^(p)).*(one_vect-u_vect))
        #return u_vect
    #end

    # Main While loop
    print(" check 4")
    global t=0
    global wavesum=sum(u_vect .> 0.5)
    global counter=0
    while wavesum<(N^2)*percentage_stop

        global t=t+1
        previous_wavesum=wavesum
        mat_pro=A*u_vect
        global u_vect=Runge_Kutta_4th_vect(mat_pro, u_vect)
        
        #global u_vect=fwd_eul_vect(A*u_vect, u_vect)
        global wavesum=sum(u_vect .> 0.5)
        global stopper=0
        m_lattice=transpose(reshape(u_vect, (N, N)))
        # Heatmap of Propagation
        #if mod(t,100)==0
            #he_25=heatmap(1:size(m_lattice,1), 1:size(m_lattice,2), m_lattice, aspect_ratio=:equal, xlims=(1, N), ylims=(1,N), c=cgrad(:imola10),title="Heatmap of Square Lattice Wave Propagation: μ="*string(lmda))
            #display(he_25)
            #gui(he_25)
        #end
        

        if (wavesum>previous_wavesum)
            print(" Countdown:"*string((N^2)*percentage_stop-previous_wavesum))
            #print(" :"*string(N-previous_wavesum))
            for i in 1:length(u_vect)

                if (u_vect[i]>Wave_Threshold && Threshold_attained_time_series[i]==0)
                    global Threshold_attained_time_series[i]=h*t
                end

            end
        end

        #if mod(t,1000)==0
            #global Threshold_attained_time_Matrix=reshape(Threshold_attained_time_series, (N, N))
            #he_t=heatmap(1:size(Threshold_attained_time_Matrix,1), 1:size(Threshold_attained_time_Matrix,2), Threshold_attained_time_Matrix, aspect_ratio=:equal, xlims=(1, N), ylims=(1,N), c=cgrad(:imola10),title="Heatmap of Square Lattice time: μ="*string(lmda))
            #display(he_t)
            #gui(he_t)
        #end


    end



    m_lattice=transpose(reshape(u_vect, (N, N)))
    # Heatmap of Propagation
    he_2=heatmap(1:size(m_lattice,1), 1:size(m_lattice,2), m_lattice, aspect_ratio=:equal, xlims=(1, N), ylims=(1,N), c=cgrad(:viridis),title="Heatmap of Square Lattice Wave Propagation: Γ="*string(lmda))
    display(he_2)
    gui(he_2)
    #savefig(he_2, "2DHeatmap_of_wavespread_lambda"*string(lmda)*".png")


    # Calculating Wavespeed
    global Threshold_attained_time_Matrix=transpose(reshape(Threshold_attained_time_series, (N, N)))
    global Speed_Matrix=zeros((N, N))

    
    # Heatmap of time
    he_t=heatmap(1:size(Threshold_attained_time_Matrix,1), 1:size(Threshold_attained_time_Matrix,2), Threshold_attained_time_Matrix, aspect_ratio=:equal, xlims=(1, N), ylims=(1,N), c=cgrad(:viridis),title="Heatmap of Square Lattice time: Γ="*string(lmda))
    display(he_t)
    gui(he_t)

    # Loop through threshold time matrix and calculate the Normal speed of each point as a corresponding entry in the speed matrix

    global non_boundary_limit=N-1
    #for i in 1:non_boundary_limit
        #for j in 1:non_boundary_limit
            # X and Y speed calculation at a point
            #if(Threshold_attained_time_Matrix[i,j]>0 && Threshold_attained_time_Matrix[i+1,j]>0 && Threshold_attained_time_Matrix[i,j+1]>0)
                #global x_speed=1/abs(Threshold_attained_time_Matrix[i+1,j]-Threshold_attained_time_Matrix[i,j])
                #global y_speed=1/abs(Threshold_attained_time_Matrix[i,j+1]-Threshold_attained_time_Matrix[i,j])
                # Calculating Norm for the point (Normal direction speed)
                #global speed_norm=sqrt(x_speed^2+y_speed^2)
                # Filling in speed matrix point
                #Speed_Matrix[i,j]=speed_norm
            #end
        #end
    #end

    for i in cent_coord_x:N
        for j in cent_coord_y:N
            if (Threshold_attained_time_Matrix[i,j]!=0)
                distance_to_central_coord=sqrt((i-cent_coord_x)^2+(j-cent_coord_y)^2)
                global Speed=distance_to_central_coord/Threshold_attained_time_Matrix[i,j]
                global Speed_Matrix[i,j]=Speed
            end
        end
    end

     # Speed Heatmap
     he_3=heatmap(1:size(Speed_Matrix,1), 1:size(Speed_Matrix,2), Speed_Matrix, aspect_ratio=:equal, xlims=(1, N), ylims=(1,N), c=cgrad(:viridis), title="Heatmap of Speeds: Γ="*string(lmda))
     display(he_3)
     gui(he_3)
     savefig(he_3, "2DHeatmap_of_Speeds_lambda"*string(lmda)*".png")

     max_value, index = findmax(Speed_Matrix)
     println(" Maximum value in Speed_Matrix is: $max_value")

    # Determining wavespeed along angle trajectories from centre (new code)

    function Angle_point(phi,cent_coord_x,cent_coord_y,epsilon)
        # creating speed vect
        global Speed_Vect=Vector{Float64}()
        global Speed_Time_Vect=Vector{Float64}()
        # phi is between 0 and pi/2

        if(deg2rad(0)==phi)
             # Calculate gradient in quadrant and constant of straight line equation
             M=tan(deg2rad(90))
             #M=tan(phi)
             #constant=cent_coord_y-M*cent_coord_x
             # Now loops through quadrant
             for i in cent_coord_x:N
                 # point is non zero
                 if(Speed_Matrix[i,1]!=0)
                     push!(Speed_Vect,Speed_Matrix[i,1])
                     push!(Speed_Time_Vect, Threshold_attained_time_Matrix[i,1])
 
                     # Calculate ideal value in SLE 
                     #sle_point_value=M*1+constant
                     # Calculate euclidean distance of ideal point from current point
                     #distance_from_ideal_to_current=sqrt((1-1)^2+(j-sle_point_value)^2)
                     # Check to see if current point is within acceptable range of the ideal point
                     #if(distance_from_ideal_to_current<=epsilon)
                         #print("  LOLOLOLOLOL ")
                         #push!(Speed_Vect,Speed_Matrix[1,j])
                         #push!(Speed_Time_Vect, Threshold_attained_time_Matrix[1,j])
                     #end
                 end
             end
        end


        if(deg2rad(0)<phi<deg2rad(90))
            # Calculate gradient in quadrant and constant of straight line equation
            M=tan(phi)
            #constant=cent_coord_y-M*cent_coord_x
            print(" Angle is:"*string(rad2deg(phi)))
            constant=-M+1
            print(" CONSTANT IS: "*string(constant))
            print(" Equation is: ")
            print(string(M)*"i+"*string(constant))
            # Now loops through quadrant
            for i in cent_coord_x:N
                for j in cent_coord_y:N
                    # point is non zero
                    if(Speed_Matrix[i,j]!=0)
                        # Calculate ideal value in SLE 
                        sle_point_value=M*i+constant
                        # Calculate euclidean distance of ideal point from current point
                        distance_from_ideal_to_current=abs(j-sle_point_value)
                        # Check to see if current point is within acceptable range of the ideal point
                        if(distance_from_ideal_to_current<=epsilon)
                            push!(Speed_Vect,Speed_Matrix[i,j])
                            push!(Speed_Time_Vect, Threshold_attained_time_Matrix[i,j])
                        end
                    end
                end
            end




        end

        # phi is between pi/2 and pi
        if(deg2rad(90)<=phi<deg2rad(180))
            # Calculate gradient in quadrant and constant of straight line equation
            M=tan(deg2rad(90)-phi)
            #M=tan(phi)
            #constant=cent_coord_y-M*cent_coord_x
            # Now loops through quadrant
            for j in cent_coord_y:N
                # point is non zero
                if(Speed_Matrix[1,j]!=0)
                    push!(Speed_Vect,Speed_Matrix[1,j])
                    push!(Speed_Time_Vect, Threshold_attained_time_Matrix[1,j])

                    # Calculate ideal value in SLE 
                    #sle_point_value=M*1+constant
                    # Calculate euclidean distance of ideal point from current point
                    #distance_from_ideal_to_current=sqrt((1-1)^2+(j-sle_point_value)^2)
                    # Check to see if current point is within acceptable range of the ideal point
                    #if(distance_from_ideal_to_current<=epsilon)
                        #print("  LOLOLOLOLOL ")
                        #push!(Speed_Vect,Speed_Matrix[1,j])
                        #push!(Speed_Time_Vect, Threshold_attained_time_Matrix[1,j])
                    #end
                end
            end
        end

        global Speed_Vect_avg=mean(Speed_Vect)
        Angle_direction_info=[Speed_Vect_avg, Speed_Time_Vect, Speed_Vect]
        #return Speed_Vect_avg
        return Angle_direction_info
    end

    
    

    # loop through different angles (in radians)
    # print(Angle_point(pi/3,cent_coord_x,cent_coord_y,1))
    # 1 degree = 0.0174533 radians, we multiply the radian component by 10 to get 10 degrees=0.174533
    global Angle_Vect=[deg2rad(0), deg2rad(4.5), deg2rad(2*4.5), deg2rad(3*4.5), deg2rad(4*4.5), deg2rad(5*4.5), deg2rad(6*4.5), deg2rad(7*4.5), deg2rad(8*4.5), deg2rad(9*4.5), deg2rad(10*4.5), deg2rad(11*4.5), deg2rad(12*4.5),deg2rad(13*4.5),deg2rad(14*4.5), deg2rad(15*4.5), deg2rad(16*4.5), deg2rad(17*4.5), deg2rad(18*4.5), deg2rad(19*4.5),deg2rad(90)]
    global Angle_to_output_Vect=Vector{Float64}() 
    global Angle_of_speed_vect=Vector{Float64}()
    global epsilon=1
    # Want to sort Speed_Time_Vect in ascending order, and sort Speed_Vect by the same way so the times are still paired with original speeds
    using DataFrames
    global average_speed_by_angle_vect=[]

    for i in 1:length(Angle_Vect)


        # Getting the avarage wavespeed in a given direction
        Angle_direction_info=Angle_point(Angle_Vect[i],cent_coord_x,cent_coord_y,epsilon)
        global average_wavespeed_output=Angle_direction_info[1]
        global Speed_Time_Vect=Angle_direction_info[2]
        global Speed_Vect=Angle_direction_info[3]
        #global Time_Vect=Angle_direction_info[4]
        push!(Angle_to_output_Vect, average_wavespeed_output)
        push!(Angle_of_speed_vect, Angle_Vect[i])

        # Want to sort Speed_Time_Vect in ascending order, and sort Speed_Vect by the same way so the times are still paired with original speeds
        df=DataFrame(a=Speed_Time_Vect,b=Speed_Vect);
        df_sorted=sort(df, :a)
        Sorted_time_vector=df_sorted[!,1]
        global Sorted_Speed_vector=df_sorted[!,2]
        # Getting the evolution in time of wavespeed in a given direction - plotting Speed_Time_Vect against Speed_Vect after using Angle_point of given angle
        display(plot(Sorted_time_vector, Sorted_Speed_vector,xlabel="Time",ylabel="Speed", mode="lines",legend=false, title="Angle="*string(rad2deg(Angle_Vect[i]))))
        #savefig(plot(Sorted_time_vector, Sorted_Speed_vector,xlabel="Time",ylabel="Speed", mode="lines",legend=false, title="Angle="*string(rad2deg(Angle_Vect[i]))), "2DSquares_Speed_evolution_by_angle:"*string(Angle_Vect[i])*".png")

        replace!(Sorted_Speed_vector, Inf=>0)
        #print(Sorted_Speed_vector)
        #print(" ")
        #print(mean(Sorted_Speed_vector))
        push!(average_speed_by_angle_vect, mean(Sorted_Speed_vector))
        push!(sorted_speed_by_angle_leading_point, Sorted_Speed_vector[length(Sorted_Speed_vector)])


        # Select p(0) and p(θ)
        if i==1
            global p_0=Sorted_Speed_vector[length(Sorted_Speed_vector)]
            push!(p_θ_vect, p_0)
        end
        if i!=1
            global p_θ=Sorted_Speed_vector[length(Sorted_Speed_vector)]
            push!(p_θ_vect, p_θ)
        end

    end


    display(plot(Angle_Vect, sorted_speed_by_angle_leading_point, linewidth=5, thickness_scaling = 1, c=:red, xlabel="θ in Radians", ylabel="Speed", legend=false))
    display(hline!([0.22580876946449277], linestyle=:dash, label=["κ="*string(0.22580876946449277)], legend=false))
    #savefig(plot(Angle_Vect, sorted_speed_by_angle_leading_point, linewidth=5, thickness_scaling = 1, c=:red, xlabel="θ in Radians", ylabel="Speed", legend=false), "2DSquare_wavespeed_by_angles"*string(lmda)*".png")

    plot(Angle_Vect, average_speed_by_angle_vect, linewidth=5, thickness_scaling = 1, c=:red, xlabel="θ in Radians", ylabel="Average Speed", legend=false)
    display(hline!([0.22580876946449277], linestyle=:dash, label=["κ="*string(0.22580876946449277)], legend=false))
    #savefig(plot(Angle_Vect, average_speed_by_angle_vect, linewidth=5, thickness_scaling = 1, c=:red, xlabel="θ in Radians", ylabel="Average Speed"), "2DSquare_Avg_wavespeed_by_angles_lined"*string(lmda)*".png")


    if z==1
        #print(" PRINTO!!")
        #p_0_vect=repeat([p_0],length(Angle_Vect))
        #p_0_div_p_θ=p_θ_vect./p_0_vect
        #display(plot(Angle_Vect[1:9], p_0_div_p_θ[1:9], xlabel="θ in Radians", ylabel="P(0)/P(θ)", marker=:circle, legend=false))
    end

    if z!=1
        #print(" printing lol")
        #p_0_vect=repeat([p_0],length(Angle_Vect))
        #p_0_div_p_θ=p_θ_vect./p_0_vect
        #display(plot!(Angle_Vect[1:9], p_0_div_p_θ[1:9], xlabel="θ in Radians", ylabel="P(0)/P(θ)", marker=:circle, legend=false))
        #savefig(plot!(Angle_Vect[1:9], p_0_div_p_θ[1:9], xlabel="θ in Radians", ylabel="P(0)/P(θ)", marker=:circle, legend=false), "2d_square_graph_c_3.png")
    end

    
end

