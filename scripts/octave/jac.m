close all
clear
clc


%%




function sigma = compute_cell_stress(Ux0, Ux1, Ux2, Ux3, Uy0, Uy1, Uy2, Uy3)

def_grad = compute_cell_def_grad (Ux0, Ux1, Ux2, Ux3, Uy0, Uy1, Uy2, Uy3);

second_invariant = .5 .* sum(def_grad.*def_grad);

viscos = yield_shear_stress/sqrt(second_invariant) + 2.*mu;

sigma = [viscos*def_grad(1), viscos*def_grad(2), viscos*def_grad(4)];

end


function def_grad = compute_cell_def_grad (Ux0, Ux1, Ux2, Ux3, Uy0, Uy1, Uy2, Uy3)

h_cell = .25*(hdof(1)+hdof(2)+hdof(3)+hdof(4));
Ux_cell = .25*(Ux0+Ux1+Ux2+Ux3); 
Uy_cell = .25*(Uy0+Uy1+Uy2+Uy3); 

grad_cell_ux = [.5 * ( (Ux3/hdof[3] - Ux2/hdof(3)) + (Ux1/hdof(2) - Ux0/hdof(1)) )/Dx, .5 * ( (hdof[2]>epsilon ? Uxdof_2/hdof[2] : 0. - hdof[0]>epsilon ? Uxdof_0/hdof[0] : 0.) + (hdof[3]>epsilon ? Uxdof_3/hdof[3] : 0. - hdof[1]>epsilon ? Uxdof_1/hdof[1] : 0.) )/Dy];
grad_cell_uy = [.5 * ( (Uy3/hdof[3] - Uy2/hdof(3)) + (Uy1/hdof(2) - Uy0/hdof(1)) )/Dx, .5 * ( (hdof[2]>epsilon ? Uydof_2/hdof[2] : 0. - hdof[0]>epsilon ? Uydof_0/hdof[0] : 0.) + (hdof[3]>epsilon ? Uydof_3/hdof[3] : 0. - hdof[1]>epsilon ? Uydof_1/hdof[1] : 0.) )/Dy];



end