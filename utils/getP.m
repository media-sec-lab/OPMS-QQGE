function [P] = getP( lambda, cost )
e = exp(1);
JointProb = e .^(-lambda*cost);
A = ones(size(cost,1),1)*sum(JointProb);
P = JointProb./A;
end

