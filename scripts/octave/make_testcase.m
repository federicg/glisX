## make_testcase.m
##
## Generate an idealised inclined-plane test case for glisX:
##
##   inputs/dem.octbin.gz      variable "dem"     - bedrock elevation raster [m]
##   inputs/mask_in.octbin.gz  variable "mask_in" - initial release indicator (0/1)
##
## The raster is Ny rows x Nx cols with spacing "res" [m]. Physical coordinates
## follow the convention used by glisX: x = (col-1)*res grows eastwards,
## y = (Ny-row)*res grows northwards (raster row 1 is the northern edge).
##
## Edit the parameters below and run:
##   octave --no-gui make_testcase.m
## from the repository root (or adjust out_dir).

out_dir = "../../inputs";

res = 1.0;      # raster spacing [m]
Nx  = 101;      # number of columns
Ny  = 101;      # number of rows

slope_deg = 20;         # bed inclination along x [deg]
release_radius = 10;    # radius of the circular material release [m]

## ---------------------------------------------------------------------------

L = res * (Nx - 1);
H = res * (Ny - 1);

col = 0:(Nx - 1);
row = 0:(Ny - 1);
[X, ROW] = meshgrid(col * res, row);
Y = (Ny - 1 - ROW) * res;

## Bedrock: plane tilted along x, highest at x = 0.
dem = tand(slope_deg) * (L - X);

## Release area: circular patch centred on the domain, flagged with 1.
mask_in = double(((X - L/2).^2 + (Y - H/2).^2) <= release_radius^2);

if (exist(out_dir, "dir") != 7)
  mkdir(out_dir);
endif

save("-binary", fullfile(out_dir, "dem.octbin.gz"), "dem");
save("-binary", fullfile(out_dir, "mask_in.octbin.gz"), "mask_in");

printf("Wrote %s (%dx%d, elevation %.2f..%.2f m)\n", ...
       fullfile(out_dir, "dem.octbin.gz"), Ny, Nx, min(dem(:)), max(dem(:)));
printf("Wrote %s (%d release cells)\n", ...
       fullfile(out_dir, "mask_in.octbin.gz"), sum(mask_in(:)));
