using MKL
using Statistics
using Plots
using LinearAlgebra
using SparseArrays

# Important Variables
N=200
lmda=10
p=1
h=0.001
Wave_Threshold=0.5
global Threshold_attained_time_series=zeros(N^2)
one_vect=ones((N^2))
# Define and initialise lattice encoding vector, define stopping element
init_lattice=spzeros(N,N)
cent_coordinate=ceil(Int, N/2)
cent_coord_x=cent_coordinate
cent_coord_y=cent_coordinate
init_radius=ceil(Int, N*0.001)
init_lattice[cent_coord_x, cent_coord_y]=1
#for i in cent_coordinate-init_radius:cent_coordinate+init_radius
    #for j in cent_coordinate-init_radius:cent_coordinate+init_radius
        #init_lattice[i,j]=1
    #end
#end

print("check 1")

# Heatmap of initial lattice
he_1=heatmap(1:size(init_lattice,1), 1:size(init_lattice,2), init_lattice, c=cgrad([:blue, :white,:red, :yellow]), title="Heatmap of Square Lattice Domain: μ="*string(lmda))
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
function Runge_Kutta_4th_vect(mat_product::SparseVector{Float64, Int64}, u_vect::SparseVector{Float64, Int64})
    k_1=lmda*(mat_product)+(u_vect.^(p)).*(one_vect-u_vect)
    k_2=lmda*((mat_product+h*(k_1/2)))+((u_vect+h*(k_1/2)).^(p)).*(one_vect-(u_vect+h*(k_1/2)))
    k_3=lmda*((mat_product+h*(k_2/2)))+((u_vect+h*(k_2/2)).^(p)).*(one_vect-(u_vect+h*(k_2/2)))
    k_4=lmda*((mat_product+h*(k_3)))+((u_vect+h*(k_3)).^(p)).*(one_vect-(u_vect+h*(k_3)))
    u_vect=u_vect+(h/6)*(k_1+2*k_2+2*k_3+k_4)
    return u_vect
end

# Main While loop


print(" check 4")
global t=0
global wavesum=sum(u_vect .> 0.5)
#print(" "*string(wavesum))
while u_vect[stopping_index]<Wave_Threshold

    #print(" check 5")
    global t=t+1
    previous_wavesum=wavesum
    
    # Do matrix multiplication
    #global mat_product=A*u_vect
    #print(" check 6")
    # Now do rk_vect method
    global u_vect=Runge_Kutta_4th_vect(A*u_vect, u_vect)
    #print(" check 7")

    # Capturing when wave passes node
    global wavesum=sum(u_vect .> 0.5)
    #print(" "*string(wavesum))
    global stopper=0
    #print(" check 8")

    if (wavesum>previous_wavesum)
        for i in 1:length(u_vect)
            
            # Calculates distance to stopping index
            if (u_vect[i]>Wave_Threshold && stopper==0)
                print(" ")
                print(abs(stopping_index-i))
                global stopper=1
            end
            
            if (u_vect[i]>Wave_Threshold && Threshold_attained_time_series[i]==0)
                global Threshold_attained_time_series[i]=h*t
            end

        end
    end


end



m_lattice=reshape(u_vect, (N, N))
he_2=heatmap(1:size(m_lattice,1), 1:size(m_lattice,2), m_lattice, c=cgrad([:blue, :white,:red, :yellow]), title="Heatmap of Square Lattice Wave Propagation: μ="*string(lmda))
display(he_2)
gui(he_2)
savefig(he_2, "2DHeatmap_of_wavespread_lambda"*string(lmda)*".png")


# Calculating Wavespeed
global Threshold_attained_time_Matrix=reshape(Threshold_attained_time_series, (N, N))
global Speed_Matrix=zeros((N, N))

# Loop through threshold time matrix and calculate the Normal speed of each point as a corresponding entry in the speed matrix

global non_boundary_limit=N-1
for i in 1:non_boundary_limit
    for j in 1:non_boundary_limit
        # X and Y speed calculation at a point
        if(Threshold_attained_time_Matrix[i,j]>0 && Threshold_attained_time_Matrix[i+1,j]>0 && Threshold_attained_time_Matrix[i,j+1]>0)
            global x_speed=1/abs(Threshold_attained_time_Matrix[i+1,j]-Threshold_attained_time_Matrix[i,j])
            global y_speed=1/abs(Threshold_attained_time_Matrix[i,j+1]-Threshold_attained_time_Matrix[i,j])
            # Calculating Norm for the point (Normal direction speed)
            global speed_norm=sqrt(x_speed^2+y_speed^2)
            # Filling in speed matrix point
            Speed_Matrix[i,j]=speed_norm
        end
    end
end

he_3=heatmap(1:size(Speed_Matrix,1), 1:size(Speed_Matrix,2), Speed_Matrix, c=cgrad([:blue, :white,:red, :yellow]), title="Heatmap of Speeds: μ="*string(lmda))
display(he_3)
gui(he_3)
savefig(he_3, "2DHeatmap_of_Speeds_lambda"*string(lmda)*".png")


# Determining wavespeed along angle trajectories from centre (new code)

function Angle_point(phi,cent_coord_x,cent_coord_y,epsilon)
    # creating speed vect
    global Speed_Vect=Vector{Float64}()
    global Speed_Time_Vect=Vector{Float64}()
    # phi is between 0 and pi/2
    if(deg2rad(0)<=phi<deg2rad(90))
        # Calculate gradient in quadrant and constant of straight line equation
        M=tan(phi)
        constant=cent_coord_y-M*cent_coord_x
        # Now loops through quadrant
        for i in cent_coord_x:N
            for j in cent_coord_y:N
        
                # point is non zero
                if(Speed_Matrix[i,j]!=0)

                    # Calculate ideal value in SLE 
                    sle_point_value=M*i+constant
                    # Calculate euclidean distance of ideal point from current point
                    distance_from_ideal_to_current=sqrt((i-i)^2+(j-sle_point_value)^2)
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
        constant=cent_coord_y-M*cent_coord_x
        # Now loops through quadrant
        for i in 1:cent_coord_x
            for j in cent_coord_y:N
        
                # point is non zero
                if(Speed_Matrix[i,j]!=0)

                    # Calculate ideal value in SLE 
                    sle_point_value=M*i+constant
                    # Calculate euclidean distance of ideal point from current point
                    distance_from_ideal_to_current=sqrt((i-i)^2+(j-sle_point_value)^2)
                    # Check to see if current point is within acceptable range of the ideal point
                    if(distance_from_ideal_to_current<=epsilon)
                        push!(Speed_Vect,Speed_Matrix[i,j])
                        push!(Speed_Time_Vect, Threshold_attained_time_Matrix[i,j])
                    end
                end
            end
        end
    end

    # phi is between pi and 3/2*pi
    if(deg2rad(180)<=phi<deg2rad(270))
     # Calculate gradient in quadrant and constant of straight line equation
        M=tan(phi)
        constant=cent_coord_y-M*cent_coord_x
        # Now loops through quadrant
        for i in 1:cent_coord_x
            for j in 1:cent_coord_y
        
                # point is non zero
                if(Speed_Matrix[i,j]!=0)
                    #print(" CHECK 5! ")
                    #print(phi)
                    # Calculate ideal value in SLE 
                    sle_point_value=M*i+constant
                    # Calculate euclidean distance of ideal point from current point
                    distance_from_ideal_to_current=sqrt((i-i)^2+(j-sle_point_value)^2)
                    # Check to see if current point is within acceptable range of the ideal point
                    if(distance_from_ideal_to_current<=epsilon)
                        #print("CHECK 6! ")
                        push!(Speed_Vect,Speed_Matrix[i,j])
                        push!(Speed_Time_Vect, Threshold_attained_time_Matrix[i,j])
                    end
                end
            end
        end
    end

    # phi is between 3/2*pi and 2pi
    if(deg2rad(270)<=phi<=deg2rad(360))
        # Calculate gradient in quadrant and constant of straight line equation
        M=tan(deg2rad(90)-phi)
        constant=cent_coord_y-M*cent_coord_x
        # Now loops through quadrant
        for i in cent_coord_x:N
            for j in 1:cent_coord_y
        
                # point is non zero
                if(Speed_Matrix[i,j]!=0)

                    # Calculate ideal value in SLE 
                    sle_point_value=M*i+constant
                    # Calculate euclidean distance of ideal point from current point
                    distance_from_ideal_to_current=sqrt((i-i)^2+(j-sle_point_value)^2)
                    # Check to see if current point is within acceptable range of the ideal point
                    if(distance_from_ideal_to_current<=epsilon)
                        push!(Speed_Vect,Speed_Matrix[i,j])
                        push!(Speed_Time_Vect, Threshold_attained_time_Matrix[i,j])
                    end
                end
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
global Angle_Vect=[deg2rad(0), deg2rad(10), deg2rad(20), deg2rad(30), deg2rad(40), deg2rad(50), deg2rad(60), deg2rad(70), deg2rad(80), deg2rad(90), deg2rad(100), deg2rad(110), deg2rad(120), deg2rad(130),  deg2rad(140), deg2rad(150), deg2rad(160), deg2rad(170), deg2rad(180), deg2rad(190),deg2rad(200), deg2rad(210), deg2rad(220), deg2rad(230), deg2rad(240), deg2rad(250), deg2rad(260), deg2rad(270), deg2rad(280), deg2rad(290), deg2rad(300), deg2rad(310), deg2rad(320), deg2rad(330), deg2rad(340), deg2rad(350), deg2rad(360)]
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
    push!(Angle_to_output_Vect, average_wavespeed_output)
    push!(Angle_of_speed_vect, Angle_Vect[i])

    # Want to sort Speed_Time_Vect in ascending order, and sort Speed_Vect by the same way so the times are still paired with original speeds
    df=DataFrame(a=Speed_Time_Vect,b=Speed_Vect);
    df_sorted=sort(df, :a)
    Sorted_time_vector=df_sorted[!,1]
    Sorted_Speed_vector=df_sorted[!,2]
    # Getting the evolution in time of wavespeed in a given direction - plotting Speed_Time_Vect against Speed_Vect after using Angle_point of given angle
    display(plot(Sorted_time_vector, Sorted_Speed_vector,xlabel="Time",ylabel="Speed", mode="lines",legend=false, title="Angle="*string(rad2deg(Angle_Vect[i]))))
    savefig(plot(Sorted_time_vector, Sorted_Speed_vector,xlabel="Time",ylabel="Speed", mode="lines",legend=false, title="Angle="*string(rad2deg(Angle_Vect[i]))), "2DSquares_Speed_evolution_by_angle:"*string(Angle_Vect[i])*".png")

    replace!(Sorted_Speed_vector, Inf=>0)
    #print(Sorted_Speed_vector)
    #print(" ")
    #print(mean(Sorted_Speed_vector))
    push!(average_speed_by_angle_vect, mean(Sorted_Speed_vector))
    

end

display(plot(Angle_Vect, average_speed_by_angle_vect, xlabel="θ in Radians", ylabel="Average speed", legend=false))
savefig(plot(Angle_Vect, average_speed_by_angle_vect, xlabel="θ in Radians", ylabel="Average speed", legend=false), "2DSquare_wavespeed_by_angles"*string(lmda)*".png")