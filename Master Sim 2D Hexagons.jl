using MKL
using Statistics
using Plots
using LinearAlgebra
using SparseArrays

# Important Variables
N=200
lmda=1
p=1
h=0.001
Wave_Threshold=0.5
global Threshold_attained_time_series=zeros(N^2)
one_vect=ones((N^2))
# Define and initialise lattice encoding vector, define stopping element
init_lattice=spzeros(N,N)
cent_coordinate=ceil(Int, N/2)
init_radius=ceil(Int, N*0.01)
for i in cent_coordinate-init_radius:cent_coordinate+init_radius
    for j in cent_coordinate-init_radius:cent_coordinate+init_radius
        init_lattice[i,j]=1
    end
end

print("check 1")

# Heatmap of initial lattice
he_1=heatmap(1:size(init_lattice,1), 1:size(init_lattice,2), init_lattice, c=cgrad([:blue, :white,:red, :yellow]), title="Heatmap of Hexagonal Lattice Domain")
display(he_1)
gui(he_1)

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
v_c=ones(N^2)*-6
A[diagind(A,0)]=v_c
# Far right diagonal
v_FR=ones(N^2-N)
v_FR[1:N].=2
v_FR[1]=3



# Middle right diagonal
c_m=ones(N)
c_m[1]=0
c_m[N]=2
v_m=repeat(c_m, outer = N)
#pop!(v_m) # gets rid of last (n-1) element
v_Rm=v_m[1:(N^2-(N-1))]
#v_m[1]=0
v_Rm[2:N].=2



# Right diagonal
c_1=ones(N)
c_1[1]=2
c_1[N]=0
v_R=repeat(c_1, outer = N)
pop!(v_R) # gets rid of last element
v_R[1]=3
# Far left diagonal
v_FL=reverse(v_FR)
# Middle left diagonal
v_Lm=reverse(v_Rm)


# Left diagonal
v_L=reverse(v_R)
# Creating true matrix A
A[diagind(A,N)]=v_FR
A[diagind(A,-N)]=v_FL
A[diagind(A,1)]=v_R
A[diagind(A,-1)]=v_L
A[diagind(A,-N+1)]=v_Lm
A[diagind(A,N-1)]=v_Rm

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
he_2=heatmap(1:size(m_lattice,1), 1:size(m_lattice,2), m_lattice, c=cgrad([:blue, :white,:red, :yellow]), title="Heatmap of Square Lattice Domain")
display(he_2)
gui(he_2)



# Calculating Wavespeed
global Threshold_attained_time_Matrix=reshape(Threshold_attained_time_series, (N, N))
global Speed_Matrix=zeros((N, N))

# Loop through threshold time matrix and calculate the Normal speed of each point as a corresponding entry in the speed matrix

global non_boundary_limit=N-1
for i in 1:non_boundary_limit
    for j in 1:non_boundary_limit
        # X and Y speed calculation at a point
        if(Threshold_attained_time_Matrix[i,j]>0 && Threshold_attained_time_Matrix[i+1,j]>0 && Threshold_attained_time_Matrix[i,j+1]>0)
            global x_time_dist=1/abs(Threshold_attained_time_Matrix[i+1,j]-Threshold_attained_time_Matrix[i,j])
            global y_time_dist=1/abs(Threshold_attained_time_Matrix[i,j+1]-Threshold_attained_time_Matrix[i,j])
            # Calculating Norm for the point (Normal direction speed)
            global time_norm=sqrt(x_time_dist^2+y_time_dist^2)
            # Filling in speed matrix point
            Speed_Matrix[i,j]=time_norm
        end
    end
end

he_3=heatmap(1:size(Speed_Matrix,1), 1:size(Speed_Matrix,2), Speed_Matrix, c=cgrad([:blue, :white,:red, :yellow]), title="Heatmap of Speeds: Lmda value="*string(lmda))
display(he_3)
gui(he_3)
#savefig(he_3, "Heatmap_of_Speeds_lambda"*string(lmda)*".png")