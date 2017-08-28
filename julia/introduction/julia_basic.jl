#= Some examples are borrowed from
 "Julia by example" By Samuel Colvin &
 Official manual of Julia (https://docs.julialang.org/en/stable/manual/)
=#

# 0. Introduction
println("Hello. ","你好! ","안녕하세요?");
println(nprocs());
# Atom is useful because... it allows efficient editing
# for e.g., multiple line editing: shift+cmd+L



# 1. Functions
# function to calculate the volume of a sphere
function sphere_vol(r);
    # julia allows Unicode names (in UTF-8 encoding)
    # so either "pi" or the symbol π can be used
    return 4/3*pi*r^3;
end;
sphere_vol(3);



# 2. Directories
# Set path using string variables
homedir()
Path="/Dropbox/CompEcon_class/julia/";
cd(string(homedir(),Path));
pwd()
cd("$(homedir())/Dropbox/CompEcon_class/julia/");
pwd()



# 3. Arrays
# arrays can be initialised and edited directly:
a1 = [1; 2; 3]
a1 = [1 2 3]
b1 = [1 2 3; 4 5 6]
c1 = [a1; b1]
c1[1,2]=5
c1

# arrays can also be generated from comprehensions:
a4= [1; 3; 5; 7; 9; 11; 13; 15]
a4 = [i for i = 1:2:15]

# matrix & operations & index
m1=zeros(3,3);
m1[1,2]=1;
m1[2,2]=2;
m1
sum(m1,1) # takes the sum over the first dim
maximum(m1,2) # find the max over the 2nd dim
m1[1:2,:] # grabs the rows from 1 to 2 and grab all the columns
m2 = m1 .+ 1 # add 1 to all elements; "." is an element-wise operation
m2



# 4. Loops and Map
# for loops:
for i in 1:5;
  print(i, ", ");
end;

# while loops:
i=1;
while i<=5;
 print(i,", ");
 i+=1; # Shorthand for i = i + 1
end;

# map works by applying a function on each value of an array
# useful if functions are expensive
@time a2=zeros(5); @time for i=1:5; a2[i]=sphere_vol(i) end; # 10 allocations
@time a3 = [sphere_vol(1),sphere_vol(2),sphere_vol(3),sphere_vol(4),sphere_vol(5)]; # 5 allocations
@time a4 = map(sphere_vol, [1,2,3,4,5]); # 7 allocations



# 5. Packages
# using pacakages saves time
Pkg.add("Calculus")
Pkg.update()
Pkg.status("Calculus")
using Calculus;
derivative(x -> f(x,y), [1.0 1.0])
cos(1)

# using optim package
#=
 Rosenbrock test
 In mathematical optimization, the Rosenbrock function is a non-convex function
 used as a performance test problem for optimization algorithms
 introduced by Howard H. Rosenbrock in 1960.
 It is also known as Rosenbrock's valley or Rosenbrock's banana function.
=#
Pkg.add("Optim")
using Optim; #detailed manual @ http://julianlsolvers.github.io/Optim.jl/stable/
y=0
rosenbrock(x,y) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;
result = optimize(x->rosenbrock(x,y), zeros(2), BFGS())
result = optimize(rosenbrock, [0.0, 0.0], NelderMead())
result = optimize(rosenbrock, [0.0, 0.0],BFGS(),
               Optim.Options(g_tol = 1e-12,iterations = 1000,store_trace = true,show_trace = true));
result
function g(storage, x);
 storage[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1];
 storage[2] = 200.0 * (x[2] - x[1]^2);
end;
result = optimize(rosenbrock, g, [0.0, 0.0], LBFGS())
@show result
function h(storage, x);
 storage[1, 1] = 2.0 - 400.0 * x[2] + 1200.0 * x[1]^2;
 storage[1, 2] = -400.0 * x[1];
 storage[2, 1] = -400.0 * x[1];
 storage[2, 2] = 200.0;
end;
result = optimize(rosenbrock, g, h, [0.0, 0.0], LBFGS())
@show result
Optim.minimizer(result);



# 6. Plots
# This section only works in Atom
using Plots;
plot([cumsum(rand(500).-0.5) cumsum(rand(500).-0.5)])



# 7. Parallel Computing
# Parallel programming in Julia is built on two primitives: remote references and remote calls.
# A remote reference is an object that can be used from any process to refer to an object stored on a particular process.
# A remote call is a request by one process to call a certain function on certain arguments on another (possibly the same) process.

# Let's add some processors ("workers")
addprocs(2);
nprocs()
#rmprocs(4)

# Note that where you have the data matters
# suppose you want to construct an array [100, 200, 300, 400, 500]
a=zeros(5);
for i = 1:5;
  a[i]=i*100;
end;
println(a);
# suppose you want do construct the same array parallelly
@everywhere a=zeros(5);
@parallel for i = 1:5;
    a[i]=i*100;
end;
println(a);
@spawnat 2 println(a);
@spawnat 3 println(a);
# Not what you want...
# This happened because iterations run on different processors and do not happen in a specified order
# Any variables used inside the parallel loop will be copied and broadcast to each processor


# you can manually call in all the different "a"'s generated in different processors
# OR use Julia's SharedArray
a = SharedArray{Int64}(zeros(5));
@parallel for i = 1:5;
    a[i]=i*100;
end;
println(a);

# Using Pmap for more expensive reduction operations
a =0; for i=1:3; a+=i^2; end;
a = @parallel (+) for i=1:3; i^2; end; # this calls the function
a = sum(pmap(x->x^2,[1,2,3])); # this calls the function x->x^2 once


# Example: If time permits... Let's find π using Monte Carlo
#= theory:
Let's imagine a circle with r=1 in a 2 by 2 square.
We know that (area of circle / area of square) = π/4.
Now, if we randomly draw points from a 2 by 2 square,
(the number of draws that falls inside the circle / the total number of draws) ≈ π/4.
Then, π ≈ 4*(n_circle/n_total)
=#

# non-parallel version
function find_n_circle(n_total);
 n_circle=0;
 for i=1:n_total;
  x,y = rand(2)
  if (x^2+y^2<=1);
   n_circle+=1;
  end;
 end;
return n_circle;
end;
@time 4*find_n_circle(10^2)/10^2
@time 4*find_n_circle(10^4)/10^4
@time 4*find_n_circle(10^8)/10^8
pi

# Now define findpi in every processor to parallelize the calculation
function find_n_circle(n_total);
 n_circle=0;
 for i=1:n_total;
  x,y = rand(2)
  if (x^2+y^2<=1);
   n_circle+=1;
  end;
 end;
return n_circle;
end;
tic()
a=@spawn find_n_circle((10^8)/2);
b=@spawn find_n_circle((10^8)/2);
toc()
4*(fetch(a)+fetch(b))/10^8

# This example demonstrates a powerful and often-used parallel programming pattern.
# Many iterations run independently(!) over several processes,
# and then their results are combined using some function.

# using parallel loop
function par_find_n_circle(n_total);
 n_circle= @parallel (+) for i=1:n_total;
  x,y = rand(2);
  x^2+y^2<=1 ? 1: 0; # "A?B:C" is like if-else syntax; if a condition A is true, do B. Otherwise, C.
 end;
 return n_circle;
end;
@time 4*par_find_n_circle(10^8)/10^8

# using pmap
# @parallel assigns iterations to multiple processes whereas
# pmap applies a function to all elements in some range
# in the previous example, @parallel assigins the for loop to different processors
# and collect the results at the end using the (+) reduction
# but if one wants each processor to calculate find_n_circle(N), then use pmap
@time ans=pmap(find_n_circle, [(10^8)/2, (10^8)/2])
4*sum(ans)/10^8
# Julia's pmap() is designed for the case where each function call does a large amount of work.
# In contrast, @parallel for can handle situations where each iteration is tiny,
# perhaps merely summing two numbers.

#End of Code
