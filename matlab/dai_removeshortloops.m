% [psi_out] = dai_removeshortloops (psi_in)
%       
%    INPUT:  psi_in     = linear cell array containing the factors
%                         (psi{i} is a structure with a Member field
%                         and a P field, like a CPTAB).
%      
%    OUTPUT: psi_out    = linear cell array containing the factors of psi_in,
%                         where factors have been merged to prevent short loops
%                         of length 4 in the factor graph (i.e. loops of type
%                         var1-factor1-var2-factor2-var1).
