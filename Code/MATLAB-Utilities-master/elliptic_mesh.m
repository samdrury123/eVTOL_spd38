function xyOut = elliptic_mesh(xy,maxit)
    %
    % Author: Andrew Knesnik
    % Date: 11/29/2014
    %
    % Description:
    %   This function generates a smooth 2d mesh by solving the elliptic
    %   mesh pde.
    %
    % Inputs:
    %   xy = mxn complex, giving the location of a initial 2d mesh
    %
    %   xyOut = mxn complex, giving the smoothed 2d mesh with the same
    %           boundary locations as the initial mesh, but with a smooth
    %           interior grid distribution
    %
    
    % Default number of iterations
    if exist('maxit','var') == 0
        maxit = 1000;
    end
    
    % prelim
    [m, n] = size(xy);
    
    % inputs
    convCrit = 1e-5;% convergence criteri
    
    px = 2:n-1; e = 3:n; w = 1:n-2;
    py = 2:m-1; n = 3:m; s = 1:m-2;
    for i = 1:maxit
        xyE  = xy(py, e); xyN  = xy(n, px);
        xyW  = xy(py, w); xyS  = xy(s, px);
        xyNE = xy(n, e ); xySE = xy(s, e );
        xyNW = xy(n, w ); xySW = xy(s, w );
        
        % calculate derivatives
        xy_xi = (xyE - xyW)/2;
        xy_et = (xyN - xyS)/2;
        xy_xiet = (xyNE - xySE - xyNW + xySW)/4;
        
        % calcualte coefficients
        a = abs(xy_et).^2;
        b = real(xy_xi).*real(xy_et) + imag(xy_xi).*imag(xy_et); 
        c = abs(xy_xi).^2;
        
        % solve for the new interior mesh verts
        xyP = (a.*(xyE + xyW) + c.*(xyN + xyS) - 2*b.*(xy_xiet))./(2*(a + c));
        
        % check for convergence
        err = max(max(abs(xyP - xy(py, px))));
        if (err < convCrit)
            fprintf('----- MESH CONVERGED IN %i ITERATIONS\n', i);
            break;
        end
        
        % update the mesh verts
        xy(py, px) = xyP;
    end
    xyOut = xy;
end