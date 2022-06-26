function [edge_Ux] = mid2edge_2D_x(U)
% MID2EDGE_2D interpolates values of U currently found at cell midpoints of
% a mxn grid to midpoints of horizontal cell edges using edge_U(:,i) = (U(:,i) + U(:,i-1))/2 with
% periodic BC U(:,N+1) = U(:,1), U(N+1,:) = U(1,:).
%
% Input Arguments:
% U0    = matrix size mxn; x = rows, y = cols
%
% Output Arguments:
% edge_Uy    = mxn matrix for cell edge values of U
%
% Example Usage
% U = [1 2 3; 5 8 9;; 4 1 1; 7 4 11]
% [edge_U] = mid2edge_2D_x(U);

[Nx,Ny] = size(U);
edge_Ux = zeros(Nx,Ny);

% For each col 1...Ny, shift values at x-midpoints to edges
for col = 1:Ny
    % Boundary condition   
    edge_Ux(1,col) = ( U(Nx,col) + U(1,col) )/2;
    
    for row = 2:Nx
        edge_Ux(row,col) = ( U(row,col) + U(row-1,col) )/2;
    end
end

