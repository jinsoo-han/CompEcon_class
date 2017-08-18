#= Some examples are borrowed from
 "Julia by example" By Samuel Colvin &
 Official manual of Julia (https://docs.julialang.org/en/stable/manual/)
=#

# 0. Introduction
println("Hello. ","你好! ","안녕하세요?");
println(nprocs());



# 1. Functions
# function to calculate the volume of a sphere
function sphere_vol(r);
    # julia allows Unicode names (in UTF-8 encoding)
    # so either "pi" or the symbol π can be used
    return 4/3*pi*r^3;
end;
sphere_vol(3);



# 2. String
# Set path using string variables
homedir()
Path="/Dropbox/CompEcon_class/julia/";
cd(string(homedir(),Path));
pwd()
cd("$(homedir())/Dropbox/CompEcon_class/julia/");
pwd()



# 3. Arrays
# arrays can be initialised and edited directly:
a1 = [1, 2, 3];
a1 = [1; 2; 3];
a1 = [1 2 3];
b1 = [1 2 3; 4 5 6]
c1 = [a1; b1];
println(a1)

# commans are for enumeration and semi-colons are for concatenation
a2 = [a1; 4]
println(a2)
size(b1)
a3 = [a1; 4];
println(a3)
size(a3)
a4=rand(2,2);
size(a4)
a5=rand(1,2,3)
size(a5,3)
whos()
a5=[]

# arrays can also be generated from comprehensions:
 #try to use comprehensions as they are super-fast
a4= [1, 3, 5, 7, 9]
a4 = [i for i = 1:2:10]'
a4 =inv(eye(2))
println(a4)
size(a4)

# matrix & matrix operations
m1=eye(2);
m1[2,1]=2; # assign a new value to a certain cell
m1

maximum(m1,1)    # find the max elt along dim 1

m1[1,:]


sum(m1,2)  # takes the sum over the second dimension (i.e., row)
m2 = m1 .+ 1;
m3 = m1 .+ [1 2]

 # add 1 to all elements; "." is an element-wise operation
m2



# 4. Loops and Map
# for loops:
for i in 1:5; print(i, ", "); end;

# while loops:
i=1;
while i<=5;
 print(i,", ");
 i+=1; # Shorthand for i = i + 1
end;

# like python enumerate can be used to get both the index and value in a loop:
a1 = ["A" "B"; "C" "D"];
for (iter, value) in enumerate(a1[1,:]);
 print(iter, ": ", value, ", "); # Note that Julia index starts at 1 not 0
end;

# map works by applying a function on each value of an array
a2=zeros(5);
@time for i=1:5; a2[i]=sphere_vol(i) end;
@time a3 = map(sphere_vol, [1,2,3,4,5]);
@time a4 = [sphere_vol(1),sphere_vol(2),sphere_vol(3),sphere_vol(4),sphere_vol(5)];
a5 = map(x -> x^2, [1,2,3,4,5]); #using map with an anoymous function
println(a5)



# 5. Packages
# using pacakages saves time
Pkg.add("Calculus")
Pkg.status("Calculus")
Pkg.update()
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
addprocs(4);
nprocs()
#rmprocs(4)

# remotecall is a function to call a specific worker
r = remotecall(rand, 2, 3, 3); # remotecall is telling worker 2 to create 3x3 random matrix
println(r)

# spawn is a function to call a worker
r = @spawn rand(1000000000000, 1000000000000); # spawn is telling a worker to create 3x3 random matrix
println(r)

# fetch sends the work back to the master
s = r.+1; # Note that this will not run because data has not been fetched back to the master
s = fetch(r).+1;

# sending data requires some planning
# method 1
A = rand(1000,1000);
A2 = @spawn A^2;
fetch(A2);
# method 2
B = @spawn rand(1000,1000)^2;
fetch(B)
# Both methods do the same task, but the method 2 moves the data once
# But, you could also imagine a case where method 1 is preferred (i.e., the master absolutely needs the matrix A)

# check where things are
C = rand(1000,1000)^2;
remotecall_fetch(()->C, 2);
@spawnat 2 whos()

# Note that where you have the data matters
# suppose you want to add 1 through 50
a=[0];
for i = 1:100;
  a[1]+=i;
end;
println(a);
# suppose you want do the addition parallelly
@everywhere a=[0];
@parallel for i = 1:100;
    a[1]+=i;
end;
@spawnat 2 println(a)
@spawnat 3 println(a)
# Not what you want...
# you need to either call in all the different "a"'s generated in different processors
# OR use Julia's SharedArray
a = SharedArray{Int64}(1);
@parallel for i = 1:100;
    a[1]+=i;
end;
a

# Let's find π using Monte Carlo
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
