function msh = make_quadtree_connectivity (x, y);

  nx = numel (x);
  ny = numel (y);
  
  [XX,YY] = meshgrid (x, y);
  msh.p = [XX(:), YY(:)].';

  iiv(ny,nx) = 0;
  iiv(:) = 1 : nx*ny;
  iiv(end,:) = [];
  iiv(:,end) = [];
  iiv = iiv(:).';

  msh.t = [iiv; iiv+ny; iiv+ny+1; iiv+1];
  msh.t(5,:) = 1;

endfunction

