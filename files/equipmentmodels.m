% Waterside equipment model parameters.
data = struct();

% Common temperature, constants, etc.
common = struct();
common.Tchws = 7; % *C
common.Tchwr = 12; % *C
common.Tcws = 29.44; % *C
common.Tcwr = 35; % *C
common.Twb = 26.67; % *C
common.rho = 1000; % kg/m^3
common.Cp = 4.184; % kJ/(kg *C)
common.Tref = common.Tchwr; % *C

% Chillers.
chiller = struct();
chiller.a = [1.18090; -254.56150; 0.000673480];
chiller.Qnom = 10551; % kW
chiller.Qmin = 0.2*chiller.Qnom; % kW
chiller.Qmax = 1.2*chiller.Qnom; % kW
chiller.Tchws = common.Tchws; % *C
chiller.Tchwr = common.Tchwr; % *C
chiller.Tcws = common.Tcws; % *C
chiller.rhoCp = common.rho*common.Cp; % kJ/(m^3 *C)

function [W, V, Qc] = chillermodel(Q, par)
    % Gordon-ng model.
    par.Tchws = par.Tchws + 273; % Convert from *C to K.
    par.Tchwr = par.Tchwr + 273;
    par.Tcws = par.Tcws + 273;
    W = (Q + par.a(1)/par.Tchws + par.a(2)*(1 - par.Tchws/par.Tcws)) ...
        .*(par.Tcws./(par.Tchws + par.a(3)*Q)) - Q;
    V = Q./(par.rhoCp*(par.Tchwr - par.Tchws));
    Qc = Q + W; % Cooling water load.
end%function

figure()
Q = linspace(chiller.Qmin, chiller.Qmax, 101);
[W, V, Qc] = chillermodel(Q, chiller);
subplot(3, 1, 1);
title('Chiller Model');
plot(Q, W, '-k');
ylabel('W (kW)');
subplot(3, 1, 2);
plot(Q, V, '-k');
ylabel('V (m^3/s)');
subplot(3, 1, 3);
plot(Q, Qc, '-k');
ylabel('Qc (kW)');
xlabel('Q (kW)');

% Pumps. Empirical model.
pump = struct();
pump.b = [19.9767073197906, 1000, -14.3598440800279, -21.4982225348288];
pump.Vnom = 0.7571;  % m^3/s
pump.Vmin = 0.2*pump.Vnom; % m^3/s
pump.Vmax = 1.2*pump.Vnom; % m^3/s

function W = pumpmodel(V, par)
    % Empirical pump model.
    V = V/par.Vnom; % Scale volume.
    W = par.b(1)*log(1 + par.b(2)*V) + par.b(3)*V + par.b(4);
end%function

figure()
V = linspace(pump.Vmin, pump.Vmax, 101);
plot(V, pumpmodel(V, pump), '-k');
xlabel('V (m^3/s)');
ylabel('W (kW)');
title('Pump Model');

% Cooling towers.
tower = struct();
tower.c = [4.43; 1.12; 1.11];
tower.Qnom = 2461.83; % kW
tower.scale = 4;
tower.Qmin = 0.5*tower.Qnom*tower.scale;
tower.Qmax = 1.25*tower.Qnom*tower.scale;
tower.kappa = 1.8666e-05; % kJ/kg^3
tower.Twb = common.Twb; % *C
tower.Tcws = common.Tcws; % *C
tower.Tcwr = common.Tcwr; % *C
tower.Cp = common.Cp; % kJ/(kg *C)

function W = towermodel(Q, par)
    % Tower model.
    Q = Q/par.scale;
    mw = Q/(par.Cp*(par.Tcwr - par.Tcws));
    ma = par.c(2)^(1/par.c(3))*(par.c(1)*(par.Tcwr - par.Twb)./Q ...
         - mw.^(-par.c(3))).^(-1/par.c(3));
    W = par.scale*par.kappa*ma.^3;
end%function

Q = linspace(tower.Qmin, tower.Qmax, 101);
figure();
plot(Q, towermodel(Q, tower), '-k');
xlabel('Q (kW)');
ylabel('W (kW)');
title('Cooling Tower Model');

% Storage tank parameters.
tank = StorageTank('Tchwr', common.Tchwr, 'Tchws', common.Tchws, ...
                   'Tref', common.Tref, 'uvar', 'enthalpy');
Vcold0 = 0.01*tank.Vmax;
Hcold0 = tank.rhocp*Vcold0*(tank.Tchws - tank.Tref);
Vhot0 = tank.Vmax - Vcold0; 
Hhot0 = tank.rhocp*Vhot0*(tank.Tchwr - tank.Tref);

x0 = [Hcold0; Vcold0; Hhot0; Vhot0];

tank = StorageTank('uvar', 'enthalpy'); % Simulate using enthalpy.
Nsim = 240;
t = (1:Nsim)*tank.Delta;

urand = 2*(rand(size(t)) - 1);
usin = 2*sin(2*pi()*t/24);
uask = 0.125*tank.Vmax*tank.rhocp*(tank.Tchws - tank.Tref)*(urand + usin);
[x, uget] = tank.simulate(x0, uask);
tank.plot(x);
subplot(3, 1, 1);
title('Storage Tank Model');

tank = struct('Vtot', tank.Vmax, 'K', tank.k); % Save key parameters.

% Save composite data.
data.common = common;
data.chiller = chiller;
data.pump = pump;
data.tower = tower;
save('-v7', 'equipmentmodels.mat', '-struct', 'data');

