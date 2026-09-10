close all
clear
clc

%%
g = 9.81;

L = 6.6;
H = 1.;

nu = 3.7;

t_c = (L/H)^2*nu/g/H;

xf1 = @(t) .248.*(t/t_c).^.5;
xf2 = @(t) 1.133*(t/t_c + 1.221).^.2-1;

xf = @(t) L*(xf1(t).*(t<(2.5*t_c)) + xf2(t).*(t>=(2.5*t_c)));

t = linspace(0,500,1000);

figure()
plot(t,xf(t),LineWidth=2)
grid on


%% TG2
clear

g  = 9.81;
L = 75;
x   = linspace (0, L, 2^9+1) .';
dx  = diff (x);

global diff_coeff

mass = [dx(1)/2;dx(2:end);dx(end)/2];


T   = 100;
cfl = .9;


h0 = 1+.1*exp(-0.5*( (x-L/2).^2 )/(0.2*L/2)^2); %1+1*exp(-0.5*( (x-L/2).^2 )/(0.2*L/2)^2); %(x<=L/2)*70 + (x>L/2)*7; %8 + sin (pi * x / 800);
h0 = (x<=L/2)*2 + (x>L/2)*1;
h0 = (x<=6.6)*1 + (x>6.6)*0;
% h0 = (x<=5)*1 + (x>5)*0;
u0 = zeros (size (x));

h = h0;
u = u0;

nu = 3.7;%3.7;%1.16e-3;

toll_h = 1e-5;

gamma_12 = zeros(length(u)-1,1);
gamma = zeros(length(u),1);

vel = zeros(size(h));


tic
figure()
time = 0;
while (time < T)
    % low order solution computation
    % 1st step,
    vel(:) = 0.;
    vel(h>toll_h) = u(h>toll_h)./h(h>toll_h);
    diff_coeff = sqrt(g*h) + abs(vel);
    dt = min(cfl*dx(1)/max(diff_coeff), T-time);
    


    F_h = fluxh(h,u); F_u = fluxu(h,u); 
    h_cell = (h(2:end)+h(1:end-1))/2.; u_cell = (u(2:end)+u(1:end-1))/2.;
    
    temp_der_h = (F_h(2:end)-F_h(1:end-1))./dx;
    temp_der_u = (F_u(2:end)-F_u(1:end-1))./dx;
    
    h_12 = h_cell - dt/2 .* temp_der_h; 

    gamma_12(:) = 0.;
    gamma_12(h_12>toll_h) = 3.*nu./h_12(h_12>toll_h); %3.*nu./h_12(h_12>toll_h); %k_s./(1+k_s.*h_12(h_12>toll_h)./3./nu); %3.*nu./h_12;
    u_12 = (u_cell - dt/2 .* temp_der_u)./(1+dt/2.*gamma_12./h_12); 
    u_12(isnan(u_12)) = 0.;
     
    % 2 step, low order sol
    uold = u;
    hold = h;

    diff_coeff = .5*(diff_coeff(2:end)+diff_coeff(1:end-1));
    dF_h = .5 * (h(2:end)-h(1:end-1)).*diff_coeff; dF_u = .5 * (u(2:end)-u(1:end-1)).*diff_coeff;
    F_hl = [fluxh(h_12(1),u_12(1));fluxh(h_12,u_12)-dF_h;fluxh(h_12(end),u_12(end))]; F_ul = [fluxu(h_12(1),u_12(1));fluxu(h_12,u_12)-dF_u;fluxu(h_12(end),u_12(end))];
    h = h - dt./mass .* (F_hl(2:end)-F_hl(1:end-1));
    u = u - dt./mass .* (F_ul(2:end)-F_ul(1:end-1));
    
    F_h = [fluxh(h_12(1),u_12(1));fluxh(h_12,u_12);fluxh(h_12(end),u_12(end))]-F_hl; F_u = [fluxu(h_12(1),u_12(1));fluxu(h_12,u_12);fluxu(h_12(end),u_12(end))]-F_ul;
%     h_lim = flux_limiter(h.*mass./dt, F_h);
%     u_lim = flux_limiter(u.*mass./dt, F_u);
    diff_coeff = sqrt(g*h) + abs(vel);
    h_lim = flux_limiter(h, F_h);
    u_lim = flux_limiter(u, F_u);
    h = h + dt./mass .* h_lim.*(-F_h(2:end) + F_h(1:end-1));
    u = u + dt./mass .* u_lim.*(-F_u(2:end) + F_u(1:end-1));
    
    h(h<0)=0;
    
    gamma(:) = 0.;
    gamma(h>toll_h) = 3.*nu./h(h>toll_h); %k_s./(1+k_s.*h./3./nu); %3.*nu./h;
    u = (u - .5*dt.*uold.*gamma./hold)./(1+dt/2.*gamma./h); 
    u(isnan(u)) = 0;
    u(1)=0; 
 
    time = time + dt;
    time
    plot(x,[h,u],'--o')
    grid on
%     axis equal
    drawnow
end
toc



%% utils




function Fhx = fluxh (h, U)
Fhx = U;

% Fhx = -u;
end

function Fux = fluxu (h, U)
g = 9.81;
toll_h = 1e-3;
vel = zeros(size(h));
vel(h>toll_h) = U(h>toll_h)./h(h>toll_h);
Fux = U.*vel + (1/2) * g * h.^2;

% Fux = -h;
end



function lim = flux_limiter (h, F)
global diff_coeff 
 
hmin = min(min(h, [h(2:end);h(end)]), [h(1);h(1:end-1)]);
hmax = max(max(h, [h(2:end);h(end)]), [h(1);h(1:end-1)]);

F_mat = [-F(2:end), F(1:end-1)];
P_p = sum(max(0,F_mat),2);
P_m = sum(min(0,F_mat),2);

Q_p = (hmax-h).*diff_coeff;
Q_m = (hmin-h).*diff_coeff;

R_p = (P_p==0).*1 + (P_p~=0).*min(1,Q_p./(P_p+eps));
R_m = (P_m==0).*1 + (P_m~=0).*min(1,Q_m./(P_m+eps));
 
lim = min(R_p.*(F_mat>=0) + R_m.*(F_mat<0),[],2);
  
% 
% hmin = min(min(hcell, [hcell(2:end);hcell(end)]), [hcell(1);hcell(1:end-1)]);
% hmax = max(max(hcell, [hcell(2:end);hcell(end)]), [hcell(1);hcell(1:end-1)]);
% phi = (abs(hrec-hcell)<=toll) + ...
%     ((hrec-hcell)<-toll) .* min(1, (hmin-hcell)./(hrec-hcell+eps)) + ...
%     ((hrec-hcell)> toll) .* min(1, (hmax-hcell)./(hrec-hcell+eps));
% 
% phi = min(phi,[],2);
% 
% 
% ratio = (hcell-[hcell(1);hcell(1:end-1)]) ./ ([hcell(2:end);hcell(end)]-hcell);
% phi = max(0, min(1,ratio));

end


function m = minmod(a,b,c)
m = (1/3.)*(sign(a)+sign(b)+sign(c)) .* min(min(abs(a),abs(b)), abs(c));
end

function m = minmod_modified_tvb(a,b,c,dx,K)
m = zeros(size(a));
 
ind = abs(a)<K*dx*dx;
m(ind) = a(ind);

ind = abs(a)>=K*dx*dx & sign(a)==sign(b) & sign(a)==sign(c);
m(ind) = minmod(a(ind),b(ind),c(ind));

end


function f = test(x)
f = x*x;
end

