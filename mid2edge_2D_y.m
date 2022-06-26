function [edge_Uy] = mid2edge_2D_y(U)
% MID2EDGE_2D interpolates values of U currently found at cell midpoints of
% a mxn grid to midpoints of vertical cell edges using edge_U(i,:) = (U(i,:) + U(i-1,:))/2 with
% periodic BC U(N+1,:) = U(1,:).
% here U(i,j) represents U_{i,i+1/2} on mesh
%
% Input Arguments:
% U0    = matrix size mxn; x = rows, y = cols
%
% Output Arguments:
% edge_U    = mxn matrix for cell edge values of U
%
% Example Usage
% U = [2 3 1;6 5 3; 5 8 6];
% [edge_U] = mid2edge_2D(U);

[Nx,Ny] = size(U);

for row = 1:Nx
    % Boundary condition   
    edge_Uy(row,1) = ( U(row,Ny) + U(row,1) )/2;
    for col= 2:Ny
        edge_Uy(row,col) = ( U(row,col) + U(row,col-1) )/2;
    end
end

