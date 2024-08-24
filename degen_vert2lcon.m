function [A, b, A_eq, b_eq] = degen_vert2lcon(Points)
%DEGEN_VERT2LCON Find a linearly constrained convex hull of a set of points.
%   In the case that the convex hull lies in a degenerate affine subspace, 
%   identify and return linear equality constraints (hyperplanes) 
%   characterizing said subspace.
%   
%   USAGE:
%   [A, b, A_eq, b_eq] = degen_vert2lcon(Points)
%
%   OUTPUT: 
%   [A_eq, b_eq] -- 
%       A_eq \in R^[n0 x d]
%       b_eq \in R^[n0 x 1]
%       Linear equality constraints that characterize the subspace on which
%       the convex hull lies. 
%       x is in the releant subspace iff A_eq*x - b_eq == 0 
%   [A, b] -- 
%       A \in R^{n0 x d} 
%       b \in R^{n0 x 1}
%       Linear inequality constraints that characterize the convex hull. 
%       x is in the convex hull iff A*x < b && A_eq*x == b_eq
%       
%   INPUTS:
%       Points -- dxN matrix where each column is a point in R^d.
    
    d = size(Points, 1); % dimensionality of the space the points are in. 
    
    x0 = Points(:,1);
    X = Points(:,2:end) - x0; 
    % By performing Singular Value Decomposition on these points we can
    % find axes in which the convex hull is degenerate. 

    [U,S,~] = svd(X);

    if(size(S,2) < d)
        % Correct for situations where there are less points than there are
        % dimensions.
        S = [S, zeros(d, d-size(S,2))];
    end

    if(all(S<=1e-6))
        % Catch the case where there is only one unique point in Points 
        A = [];
        b = [];
        A_eq = eye(d);
        b_eq = x0;
        return;
    end
    
    s = diag(S); % Get the singular values.


    % Columns of U corresponding to non-zero singular values form the basis
    % of the affine subspace. Columns corresponding to 0 values are
    % orthogonal to the convex hull, and can be used to form equality
    % constraints.
    convhull_basis = U(:, ~(s<=1e-6));
    orth_to_convhull = U(:, (s<=1e-6));
    
    % if x is in the convex hull, then A_eq*(x -x0) = 0 
    A_eq = orth_to_convhull';
    b_eq = A_eq * x0;   

    % for every p in the convex hull of Points, 
    % there exists l s.t p=convhull_basis*l + x0. 
    % The projection transformation gets l from p.
    convhull_projection = pinv(convhull_basis);
    L = convhull_projection * (Points - x0);
    
    % Al*l-bl <= 0 defines the convex hull of the reduced representation of
    % points (L) 
    [Al, bl] = vert2lcon(L');

    % x is _not_ in the convex hull if (convhull_basis*Al')'*(x-x0) > bl
    A = Al*convhull_basis';
    b = bl + A*x0; 
    
    % When a point (x) is in the affine subspace of the convex hull 
    % (A_eq*x-b_eq == 0), then {x | (convhull_basis*Al')'*(x-x0) - bl <= 0} 
    % captures the convex hull 
end

