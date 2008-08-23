% [logZ,q,md] = dai (psi, method, opts)
%
%    INPUT:  psi        = linear cell array containing the factors
%                         psi{i} should be a structure with a Member field
%                         and a P field, like a CPTAB).
%            method     = name of the method (see README)
%            opts       = string of options (see README)
%
%    OUTPUT: logZ       = approximation of the logarithm of the partition sum.
%            q          = linear cell array containing all final beliefs.
%            md         = maxdiff (final linf-dist between new and old single node beliefs).
