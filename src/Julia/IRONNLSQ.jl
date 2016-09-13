# ==============================================================================
# IRONNLSQ.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Performs the reweighted IRONNLSQ of the grant
# ==============================================================================
include("lsqnonneg.jl")

function IRONNLSQ(Aaux,s,q,nIter,myeps,I)
#	q = 1/100;  # diversity index
#	nIter = 3;  # Number of iterations to do for
#	myeps = 0.00001;  # reweighting epsilon?!
#	I = eye(N);

	#x = lsqnonneg(Aaux,s);
	#return x

# Removed for debug
#	A = Aaux[setdiff(1:end, 1), :];  # delete first row of 1's
	N = size(A,2);  # Number of columns
	weightI = ones(1,N);  # initial weights are 1
	for dummy=1:2
		#Atilde = [weightI; lambda*A];  # Stick the weights back in
		Aaux[1,:] = weightI;
		#xI = lsqnonneg(Atilde,s);
		x = lsqnonneg(Aaux,s);
		weightI = (I*x+myeps)'.^(q-1);  # Do the reweighting
		print("inner loop\n")
		print("$(sum(x))\n")
	end
	print("out loop\n")
	return x
end