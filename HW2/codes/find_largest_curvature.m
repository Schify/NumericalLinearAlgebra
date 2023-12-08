function [ind, lambda]=find_largest_curvature(x, y, lambdas)
    %finds the largest curvature point where the normal has a positive x
    %component.
    inds = (1:length(x))';
    Vertices = [log(x),log(y)];
    k = LineCurvature2D(Vertices);
    [~, rel_ind] = max(k);
    ind = inds(rel_ind);
    lambda = lambdas(ind);
end