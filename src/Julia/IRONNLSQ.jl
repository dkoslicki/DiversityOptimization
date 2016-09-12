# ==============================================================================
# IRONNLSQ.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Performs the reweighted IRONNLSQ of the grant
# ==============================================================================
include("lsqnonneg.jl")

function IRONNLSQ(Aaux,s)
	q = 1/100;  # diversity index
	nIter = 3;  # Number of iterations to do for
	myeps = 0.00001;  # reweighting epsilon?!
	N = size(A,2);  # Number of columns
	I = eye(N);
	A = Aaux[setdiff(1:end, 1), :];  # delete first row of 1's
	weightI = ones(1,N);  # Z=I
	for n=1:nIter
		Atilde = [weightI; lambda*A];  # Stick the weights back in
		xI = lsqnonneg(Atilde,s);
		weightI = (I*xI+myeps)'.^(q-1);  # Do the reweighting
	end
	return xI
end