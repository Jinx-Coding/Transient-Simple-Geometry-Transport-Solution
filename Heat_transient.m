% Transient heat transfer analytic solution %
% Giovanni Correra - 02/2024 %

clc
close all
clear variables

OPTIONS = optimset('Display','off','MaxIter',1e20,'MaxFunEvals',1e20...
    ,'TolFun',1e-10,'Algorithm','levenberg-marquardt');

% ------------------------------- Data --------------------------------- %

% Average values option %

avg = false;

% Penetration layer %

penet = true;

% Mass transport %

mt = false;

% Geometry %
% Slab : geometry = 1 %
% Cylinder : geometry = 2 %
% Sphere : geometry = 3 %

geometry = 1;

if geometry == 1
    fprintf('Geometry = Slab\n')
elseif geometry == 2
    fprintf('Geometry = Cylindrical\n')
elseif geometry == 3
    fprintf('Geometry = Spherical\n')
end

fprintf('\n')

% Specific length %
% Slab oh 2L height %
% Cylinder of L radius (generally valid if h > r) %
% Sphere of L radius %

L = 5/200; % (m) %

% Simulation time %

tsim = 60; % (s) 86400 days%
dt = 10; % (s) %

% Spherical approximation %

approx = false;

m = 0.3; % (kg) % 

if approx == true
    V = m/rho;
    L = 0.5 * (6*V/pi)^(1/3);
end

% Generic data %

rho = 500; %(kg/m3) %

% Thermal data %

T_inf = 250 + 273.15; % (K) %
T0 = 25 + 273.15; % (K) %
k = 0.429; % (W/(mK)) %      
h = 30; % (W/(m2K)), for boiling water  h = 1e4 (W/(m2K)) %
cp = 2210; % (J/(kgK) %

% Mass transfer data %

om_inf = 0.2475; % (kg/kg) %
om0 = 0.00052; % (kg/kg) %
Keq = 1.57e-5; % (-), interface equilibrium constant %
kc = 0.01; % (m/s) %
Diff = 5e-10; % (m2/s) %

% Position (x = 0 means center) %

x = 2/1000; % (m) %

if x == 0
    x = 1e-20;
end

% ------------------------ General solution ---------------------------- %

alpha = therdiff(k,rho,cp);

% Time step evaluation %

t = linspace(0,tsim,round(tsim/dt)+1);

Fo = zeros(1,length(t));
expn = zeros(1,length(t));
E = zeros(1,length(t));
teta0 = zeros(1,length(t));
teta = zeros(1,length(t));
d = zeros(1,length(t));
tetat = zeros(1,length(t));
ft = zeros(1,length(t));

for i = 1 : length(t)
    
    d(i) = penetration(t(i),alpha);

    if d(i) <= L && geometry == 1 && penet == true
        Bi = biot(h,d(i)/2,k,Keq,kc,Diff,mt);
    else
        Bi = biot(h,L,k,Keq,kc,Diff,mt);
    end

    [lambda1,A1] = eigenvalues(Bi,geometry);

    Fo(i) = fourier(t(i),alpha,L,Diff,mt);
    expn(i) = exponential(A1,lambda1,Fo(i));
    E(i) = error(A1,Bi,lambda1,Fo(i),geometry);
    teta0(i) = theta0(expn(i),E(i));
    teta(i) = theta(teta0(i),lambda1,x,L,geometry);
    tetat(i) = tetatotal(A1,lambda1,Fo(i),Bi,x,L,geometry);
    if mt == true
        ft(i) = final(tetat(i),om0,om_inf);
    else
        ft(i) = final(tetat(i),T0,T_inf) - 273.15;
    end

end

% ---------------------- Average solution ------------------------------ %

tetaw = zeros(1,length(t));
tetac = zeros(1,length(t));
tetah = zeros(1,length(t));
Th = zeros(1,length(t));
Tc = zeros(1,length(t));
Tw = zeros(1,length(t));
Tavg = zeros(1,length(t));

tetafc = theta(teta0(1),lambda1,0,L,geometry); % First value at center %
Tfc = final(tetafc,T0,T_inf);
tetafw = theta(teta0(1),lambda1,L,L,geometry); % First value at wall %
Tfw = final(tetafw,T0,T_inf);
Tf = (Tfw+Tfc)/2; % First guess value %

if avg == true
    for i = 1 : length(t)
        tetaw(i) = theta(teta0(i),lambda1,L,L,geometry);
        Tw(i) = final(tetaw(i),T0,T_inf);
        tetac(i) = theta(teta0(i),lambda1,0,L,geometry);
        Tc(i) = final(tetac(i),T0,T_inf);
        tetah(i) = theta(teta0(i),lambda1,L/2,L,geometry);
        Th(i) = final(tetah(i),T0,T_inf);
        n = fsolve(@(n) average(alpha,L,Tw(i),T0,t(i),n,geometry), ...
            Tf,OPTIONS);
        if n > Tw(i)
            n = Tw(i);
        elseif n < Tc(i)
            n = Tc(i);
        end
        Tavg(i) = n;
        Tf = n;
    end
end

% ----------------------- Post - Processing ---------------------------- %

fprintf('Bi = %.3f (-),',Bi)
fprintf('   lambda1 = %.3f (-),',lambda1)
fprintf('   A1 = %.3f (-)\n',A1)
fprintf('\n')
fprintf('x = %.4f (m),',x)
fprintf('   position = %.2f (-),',x/L)
fprintf('   alpha = %.4f * 1e-6 (m2/s)\n',alpha*1e6)
fprintf('\n')

for i = 1 : length(t)

    if t(length(t))<600
        fprintf('t = %.2f (s)   ',t(i))
    elseif t(length(t))<3600 && t(length(t))>=600
        fprintf('t = %.2f (min)   ',t(i)/60)
    elseif t(length(t))<3*86400 && t(length(t))>=3600
        fprintf('t = %.2f (h)   ',t(i)/3600)
    elseif t(length(t))>=3*86400
        fprintf('t = %.2f (d)   ',t(i)/86400)
    end
    fprintf('Fo = %.4f (-)   ',Fo(i))
    fprintf('E = %.3f (-)   ',E(i))
    fprintf('theta0 = %.3f (-)   ',teta0(i))
    fprintf('theta = %.3f (-)   ',teta(i))
    if mt == true
        fprintf('om = %.3f %% (-)   ',ft(i)*1e2)
    else
        fprintf('T = %.3f (C)   ',ft(i))
    end

    if avg == true
        fprintf('<T> = %.3f (C)',Tavg(i)-273.15)
    end

    fprintf('\n')

end

fprintf('\n')

% ------------------ Graphical post - processing ----------------------- %

figure (1)

if geometry == 1
    sgtitle('Slab geomety')
elseif geometry == 2
    sgtitle('Cylindrical geometry')
elseif geometry == 3
    sgtitle('Spherical geometry')
end

if avg == true
    subplot(131)
else
    subplot(121)
end

plot(t,ft,'linewidth',1.5)
if mt == true
    title('Concentration')
    ylabel('om (-)')
else
    title('Temperature')
    ylabel('T (C)')
end
xlabel('t (s)')
axis([0,t(length(t)),0,ft(length(ft))])
axis square
grid on
box on

if avg == true
    subplot(132)
else
    subplot(122)
end

plot(t,E,'linewidth',1.5)
title('Error')
xlabel('t (s)')
ylabel('error (-)')
axis([0,t(length(t)),0,max(E)])
axis square
grid on
box on

if avg == true

    subplot(133)
    plot(t,Tavg-273.15,'b','linewidth',1.5)
    hold on
    plot(t,Tw-273.15,'r--','linewidth',1.2)
    hold on
    plot(t,Tc-273.15,'g--','linewidth',1.2)
    hold on
    plot(t,Th-273.15,'m--','linewidth',1.2)
    hold on
    yline(T_inf-273.15,'k--','linewidth',1.1)
    title('Average temperature vs limit temperatures')
    legend('<T>','T wall','T center','T half','T inf')
    xlabel('t (s)')
    ylabel('T (C)')
    axis square
    grid on
    box on

end

% --------------------------- Functions -------------------------------- %

% Biot Number %

function Bi = biot(h,L,k,Keq,kc,Diff,mt)

Bi = h*L/k;

if mt == true

    Bi = Keq*kc*L/Diff;

end

end

% Thermal diffusivity %

function alpha = therdiff(k,rho,cp)

alpha = k/(rho*cp);

end

% Eigenvalues %

function [lambda1,A1] = eigenvalues(Bi,geometry)

if geometry == 1
    n = 2.139;
    lambda0 = Bi^0.5;
    lambdainf = pi/2;
elseif geometry == 2
    n = 2.238;
    lambda0 = (2*Bi)^0.5;
    lambdainf = 2.4048255;
elseif geometry == 3
    n = 2.314;
    lambda0 = (3*Bi)^0.5;
    lambdainf = pi;
end

lambda1 = lambdainf/(1+(lambdainf/lambda0)^n)^(1/n);

if lambda1 == 0
    lambda1 = 1e-20;
end

if geometry == 1
    A1 = 2*sin(lambda1)/(lambda1+sin(lambda1)*cos(lambda1));
elseif geometry == 2
    J0 = besselj(0,lambda1);
    J1 = besselj(1,lambda1);
    A1 = 2*J1*lambda1/((lambda1^2)*(J0^2 + J1^2));
elseif geometry == 3
    A1 = (2*Bi*(lambda1^2 + (Bi-1))^0.5)/(lambda1^2 + Bi^2 - Bi);
end

end

% Penetration front %

function d = penetration(t,alpha)

d = 3.66*(alpha*t)^0.5;

end

% Fourier Number %

function Fo = fourier(t,alpha,L,Diff,mt)

Fo = alpha * t / L^2;

if mt == true

    Fo = Diff*t/L^2;

end

end

% Exponential part %

function expn = exponential(A1,lambda1,Fo)

expn = A1*exp(-Fo*lambda1^2);

end

% Approximation error % 

function E = error(A1,Bi,lambda1,Fo,geometry)

if geometry == 1
    c = 11;
elseif geometry == 2
    c = 15;
elseif geometry == 3
    c = 19;
end

E = ((A1-1)*Bi^(-lambda1*Fo))*exp(-c*Fo);

end

% Theta 0 %

function teta0 = theta0(expn,E)

teta0 = min(expn-E,1); % Correction if theta > 1 there is heat generation %

end

% Theta function %

function teta = theta(teta0,lambda1,x,L,geometry)

    if geometry == 1
        teta = teta0*cos(lambda1*x/L);
    elseif geometry == 2
        J0 = besselj(0,lambda1*x/L);
        teta = teta0*J0;
    elseif geometry == 3
        teta = teta0*(sin(lambda1*x/L))/(lambda1*x/L);
    end

end

% Final temperature %

function f = final(teta,f0,f_inf)

f = teta*(f0-f_inf) + f_inf;

end

% Averages zero function %

function n = average(alpha,L,Tw,T0,t,Tavg,geometry)

if geometry == 1
    n = t - (((2*L/pi)^2)/alpha)*log(8*(Tw-T0)/(Tw-Tavg));
elseif geometry == 2
    n = t - ((L^2)/(5.78*alpha))*log(0.692*(Tw-T0)/(Tw-Tavg));
elseif geometry == 3
    n = t - ((L^2)/(9.87*alpha))*log(0.608*(Tw-T0)/(Tw-Tavg));
end

end

% Overall teta function %

function tetat = tetatotal(A1,lambda1,Fo,Bi,x,L,geometry)

if geometry == 1
    c = 11;
    g = cos(lambda1*x/L);
elseif geometry == 2
    c = 15;
    g = besselj(0,lambda1*x/L);
elseif geometry == 3
    c = 19;
    g = (sin(lambda1*x/L))/(lambda1*x/L);
end

pre = A1*exp(-Fo*lambda1^2) - ((A1-1)*Bi^(-lambda1*Fo))*exp(-c*Fo);
tetat = pre * g;
tetat = min(tetat,1);

end








