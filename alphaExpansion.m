function [labels, energy, time] = alphaExpansion(unary, vertC, horC)
	[K, N, M] = size(unary);
	h = toAlphaExpansionFormat(unary, vertC, horC);
	t = cputime;
	GCO_Expansion(h);
	labels = GCO_GetLabeling(h);
	labels = reshape(labels, N, M);
	[energy, ~, ~] = GCO_ComputeEnergy(h);
	GCO_Delete(h);
	time = cputime - t;
end
