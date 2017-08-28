# Main Code: BLP Estimation;

# 0. Introduction;
	# Julia version check;
		#VERSION
	# Number of Processors;
		#addprocs();
	# Directories;
		#@everywhere UserPath = "C:/Users/jinshan";
		#@everywhere UserPath = "/Users/jinsoohan";
		#@everywhere UserPath = homedir()
		#@everywhere Path = "/Dropbox/CompEcon_class/julia";
		#@everywhere cd(string(UserPath,Path));
		#@everywhere cd("./Dropbox/CompEcon_class/julia");
		#pwd()
		#@everywhere cd(string(UserPath,Path));
	# Load standard packages;
		@everywhere using DataFrames, Optim, Distributions;
	# Load user-written packages;
		@everywhere include("$(homedir())/Dropbox/CompEcon_class/julia/blp/blp_func.jl");
		#using BLP_module
	#include(string(UserPath,Path,"/buildquadrature_v1.jl"));
		#@everywhere include(string(UserPath,Path,"/buildconsumertype_v1.jl"));

# 1. Prep Data
	# Sample size
	println("////////////////////");
	println("Specification");
	@everywhere  nbRC=2;
	println("nb RC = ",nbRC);
	println("////////////////////");

	# Load parameter vector
	@everywhere vNLparam_true=Array(readtable("$(homedir())/Dropbox/CompEcon_class/julia/blp/blp_data/blp_example_param.csv"));
	println("vNLparam_true : ",vNLparam_true);
	@everywhere vLParam=zeros(4,1);

	# Load characteristics data
	@everywhere mData=readtable("$(homedir())/Dropbox/CompEcon_class/julia/blp/blp_data/blp_example_data.csv");
		# Check data with eye
		#head(mData); #tail(mData);
	@everywhere mX=Array(mData[:,[:intercept, :x1, :x2, :x3]]);
	@everywhere vDelta_true=Array(mData[:,:delta_true]);
	@everywhere vDelta0=Array(vDelta_true);
	@everywhere mDiffIV=Array(mData[:,[:z1,:z2]]);
	@everywhere vShare=Array(mData[:,:share]);
	@everywhere vShare0=Array(mData[:,:share0]);
	@everywhere mZ=Array(mData[:,[:x2,:x3]]);
	@everywhere mIV=Array([mX mDiffIV]);
	@everywhere T=size(unique(mData[:marketid]),1);
	@everywhere J=size(unique(mData[:productid]),1);
	@everywhere n=J*T;

	# Load integration grid/weight
	@everywhere mGrid=readtable("$(homedir())/Dropbox/CompEcon_class/julia/blp/blp_data/blp_example_grid.csv");
  	@everywhere mEta=Array(mGrid[:,[:eta1,:eta2]]);
  	@everywhere vWeight=Array(mGrid[:weight]);
  	@everywhere nbGrid=size(vWeight,1);

  # Panel Structure
  @everywhere mPanelID=mData[:,[:productid,:marketid]];
  @everywhere  id=[Array(mPanelID[:marketid]).==t for t in unique(Array(mPanelID[:marketid]))];

  # Weighting matrix
  @everywhere A=inv(mIV'mIV);

	#= Test that the inversion algorithm works
	vDelta0=Array(vDelta_true + rand(Uniform(0,1),n));
	println("norm between true and starting delta: ",norm(vDelta0-vDelta_true));
	#println("inverted delta: ",vDelta0[1:5]);
	#println("true delta: ",vDelta_true[1:5]);
		#check inputs of inverse fn
		#methods(inverse);
	temp=map(t->inverse(vDelta0,vNLparam_true,t),[i for i=1:T]);
	vDelta0=zeros(n); for t =1:T; vDelta0[id[t],:]=temp[t]; end;
	println("norm between true and inverted delta: ",norm(vDelta0-vDelta_true));
	#println("inverted delta: ",vDelta0[1:5]);
	#println("true delta: ",vDelta_true[1:5]);
	=#

# 2. GMM Estimation
  # Estimate the model (non-linear parameters) using "true" starting value
	@everywhere vNLparam_est=vNLparam_true;
	@everywhere vDelta0=vDelta_true;
	res=optimize(x->gmm_obj(vDelta0,x,A), [vNLparam_est[1], vNLparam_est[2]], BFGS(), Optim.Options(show_trace=true))
	println("Convergence status: ",Optim.converged(res))
	if Optim.converged(res)==false;
	res=optimize(x->gmm_obj(vDelta0,x,A), [vNLparam_est[1], vNLparam_est[2]], NelderMead(), Optim.Options(show_trace=true))
	end;
	vNLparam_est=Optim.minimizer(res)

	# Estimate the linear parameters
	@time temp=map(t->inverse(vDelta0,vNLparam_est,t),[i for i=1:T]);
	@time temp=pmap(t->inverse(vDelta0,vNLparam_est,t),[i for i=1:T]);
	vDelta0=zeros(n); for t =1:T; vDelta0[id[t],:]=temp[t]; end;
	vLparam_est=ivreg(vDelta0,mX,mIV,A);
	vParam_full=[vNLparam_est; vLparam_est];

	# Report
	println("non-linear parameters : ",vNLparam_est);
	println("linear parameters : ",vLparam_est);

# End of the Code;
