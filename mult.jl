#=
 * Author: Thomas Johnson, thomasjohnso2020@my.fit.edu
 * Course: CSE 4250, Fall 2023
 * Project: Project #1, Multilateration
 * Implementation: 1.9.3
 =#

 using Printf


# Define the Satellite and Location structs
struct Satellite
   x_pos::Float64
   y_pos::Float64
   z_pos::Float64
end

struct Location
   time_1::Float64
   time_2::Float64
   time_3::Float64
   time_4::Float64
end

global satellites = []
global locations = []

# Misc Debug functions to make sure my input gathering
# was correct
function print_satellites(s)
   println("X:$(s.x_pos) Y:$(s.y_pos) Z:$(s.z_pos)")
end

function print_pings(l)
   println("t1:$(l.time_1) t2:$(l.time_2) t3:$(l.time_3) t4:$(l.time_4)")
end

function print_info()
   for e in satellites
      print_satellites(e)
   end
   
   for e in locations
      print_pings(e)
   end
end


# Get input and store into arrays satellites, and locations
function get_input()
   # read satellite positions
   for _ = 1:4
      line = readline()  # Read a line of input
      x, y, z = split(line)  # Split the line into x, y, and z coordinates
      push!(satellites, Satellite(parse(Float64, x), parse(Float64, y), parse(Float64, z)))
   end

   # Read time for location to ping the stallites
   while !eof(stdin)
      line = readline()  # Read a line of input
      t1, t2, t3, t4 = split(line)  # Split the line into time values
      push!(locations, Location(parse(Float64, t1), parse(Float64, t2), parse(Float64, t3), parse(Float64, t4)))
   end
end

# Function to grab sign from N for the equation of Q
function sgn(x)
   if x > 0
      return 1
   elseif x < 0
      return -1
   else
      return 0
   end
end

# All the math
function multilateration()
   # Set satellite values
   x_i::Float64, y_i::Float64, z_i::Float64 = satellites[1].x_pos, satellites[1].y_pos, satellites[1].z_pos
   x_j::Float64, y_j::Float64, z_j::Float64 = satellites[2].x_pos, satellites[2].y_pos, satellites[2].z_pos
   x_k::Float64, y_k::Float64, z_k::Float64 = satellites[3].x_pos, satellites[3].y_pos, satellites[3].z_pos
   x_l::Float64, y_l::Float64, z_l::Float64 = satellites[4].x_pos, satellites[4].y_pos, satellites[4].z_pos

   # Distances
   
   x_ij::Float64 = x_i - x_j
   x_ji::Float64 = x_j - x_i
   x_ki::Float64 = x_k - x_i
   x_jk::Float64 = x_j - x_k
   x_lk::Float64 = x_l - x_k

   y_jk::Float64 = y_j - y_k
   y_ji::Float64 = y_j - y_i
   y_ki::Float64 = y_k - y_i
   y_lk::Float64 = y_l - y_k

   z_ji::Float64 = z_j - z_i
   z_ki::Float64 = z_k - z_i
   z_jk::Float64 = z_j - z_k
   z_lk::Float64 = z_l - z_k

   # Speed of light

   c = 299792458.0

   for l in locations
      R_i::Float64 = l.time_1
      R_j::Float64 = l.time_2
      R_k::Float64 = l.time_3
      R_l::Float64 = l.time_4 

      # Convert arrival times from nanoseconds to seconds
      R_i_sec::Float64 = R_i * 1e-9
      R_j_sec::Float64 = R_j * 1e-9
      R_k_sec::Float64 = R_k * 1e-9
      R_l_sec::Float64 = R_l * 1e-9

      # Calculate the differences in seconds and then convert to meters
      R_ij::Float64 = c * (R_i_sec - R_j_sec) 
      R_ik::Float64 = c * (R_i_sec - R_k_sec)
      R_kj::Float64 = c * (R_k_sec - R_j_sec) 
      R_kl::Float64 = c * (R_k_sec - R_l_sec)

      X_ijy::Float64 = (R_ij * y_ki) - (R_ik * y_ji)
      X_ikx::Float64 = (R_ik * x_ji) - (R_ij * x_ki)
      X_ikz::Float64 = (R_ik * z_ji) - (R_ij * z_ki)
      X_kjy::Float64 = (R_kj * y_lk) - (R_kl * y_jk)
      X_klx::Float64 = (R_kl * x_jk) - (R_kj * x_lk)
      X_klz::Float64 = (R_kl * z_jk) - (R_kj * z_lk)

      Si2::Float64 = (x_i^2) + (y_i^2) + (z_i^2)
      Sj2::Float64 = (x_j^2) + (y_j^2) + (z_j^2)
      Sk2::Float64 = (x_k^2) + (y_k^2) + (z_k^2)
      Sl2::Float64 = (x_l^2) + (y_l^2) + (z_l^2)

      Rij2xyz::Float64 = (R_ij^2) + Si2 - Sj2
      Rik2xyz::Float64 = (R_ik^2) + Si2 - Sk2
      Rkj2xyz::Float64 = (R_kj^2) + Sk2 - Sj2
      Rkl2xyz::Float64 = (R_kl^2) + Sk2 - Sl2

      # Equations

      # Def ABC
      A::Float64 = X_ikx / X_ijy
      B::Float64 = X_ikz / X_ijy
      C::Float64 = X_klx / X_kjy
      D::Float64 = X_klz / X_kjy


      # Def EF
      E::Float64 = ((R_ik * Rij2xyz) - (R_ij * Rik2xyz)) / (2 * X_ijy)
      F::Float64 = ((R_kl * Rkj2xyz) - (R_kj * Rkl2xyz)) / (2 * X_kjy)

      # Def GH, IJ 

      G::Float64 = (D - B) / (A - C)
      H::Float64 = (F - E) / (A - C)

      I::Float64 = A * G + B
      J::Float64 = A * H + E

      # Def K L

      K::Float64 = Rik2xyz + (2 * x_ki * H) + (2 * y_ki * J)
      L::Float64 = 2 * ((x_ki * G) + (y_ki * I) + z_ki)

      # Def M N O

      squared_sum = G^2 + I^2 + 1
      L_squared = L^2
      K_squared = K^2
      R_ik_squared = R_ik ^ 2

      M::Float64 = ((4 * R_ik_squared) * squared_sum) - L_squared
      N::Float64 = (8 * R_ik_squared * (G * (x_i - H) + I * (y_i - J) + z_i)) + (2 * L * K)
      O::Float64 = (4 * R_ik_squared * ((x_i - H)^2 + (y_i - J)^2 + z_i^2)) - K_squared

      Q_1::Float64 = N + (sgn(N) * sqrt((N^2) - (4 * M * O)))
      Q_2::Float64 = N - (sgn(N) * sqrt((N^2) - (4 * M * O)))

      # Calculate the two possible locations

      z_1 = Q_1 / (2 * M)
      x_1 = (G*z_1) + H
      y_1 = (I*z_1) + J
      r_1 = sqrt(x_1^2 + y_1^2 + z_1^2)

      z_2 = Q_2 / (2 * M)
      x_2 = (G*z_2) + H
      y_2 = (I*z_2) + J
      r_2 = sqrt(x_2^2 + y_2^2 + z_2^2)

      # Round to nearest meter

      x_1int::Int = round(Int, x_1)
      y_1int::Int = round(Int, y_1)
      z_1int::Int = round(Int, z_1)
      r_1int::Int = round(Int, r_1)

      x_2int::Int = round(Int, x_2)
      y_2int::Int = round(Int, y_2)
      z_2int::Int = round(Int, z_2)
      r_2int::Int = round(Int, r_2)


      # Format output 

      formatted_output = "g= $(@sprintf("%4.2e", G)), h= $(@sprintf("%4.2e", H)), j= $(@sprintf("%4.2e", J)), m= $(@sprintf("%4.2e", M)), o= $(@sprintf("%4.2e\n", O))"
      print(formatted_output)


      formatted_output = "+) x= $(@sprintf("%10d", x_1int)), y= $(@sprintf("%10d", y_1int)), z= $(@sprintf("%10d", z_1int)); r= $(@sprintf("%10d", r_1int))"
      println(formatted_output)

      formatted_output = "-) x= $(@sprintf("%10d", x_2int)), y= $(@sprintf("%10d", y_2int)), z= $(@sprintf("%10d", z_2int)); r= $(@sprintf("%10d", r_2int))"
      println(formatted_output)
      println()
   end

end

function main()
   get_input()
   multilateration()
end



main()
