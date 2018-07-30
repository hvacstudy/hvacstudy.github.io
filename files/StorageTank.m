classdef StorageTank
% Container class for a nonlinear storage tank.

properties (GetAccess=public, SetAccess=protected)
    Vmax;
    rhocp;
    k;
    Tchwr;
    Tchws;
    Tref;
    Delta;
    uvar;
    safe;
    ix = struct('Hcold', 1, 'Vcold', 2, 'Hhot', 3, 'Vhot', 4); % x row names.
    iu = struct('add', 1, 'withdraw', 2); % u row names.
    integrator;
    Nx = 4;
    Nu = 2;
    safetylam = -250;
    safetyfrac = 0.001;
end%properties

methods (Access=public)
    function self = StorageTank(varargin)
        % self = StorageTank(...)
        %
        % Initialize the object with the given parameters:
        %
        % - Vmax: Tank volume in m^3
        % - rhocp: Water heat capacity in kJ/(m^3*K)
        % - k: Tank heat transfer coefficient (m^3)
        % - Tchwr, Tchws: Temperatures for chilled water return, supply (*C)
        % - Tref: Enthalpy reference temperature (*C, defaults to Tchwr)
        % - Delta: Timestep (h)
        % - uvar: 'enthalpy' or 'volume' to choose variable for u.
        % - safe: True or False whether to protect against over/underflowing
        %         the tank (default true)
        %
        % Arguments should be passed as 'key', value pairs.
        persistent parser
        if isempty(parser)
            parser = mpctools.ArgumentParser();
            parser.add('Vmax', 7.5, {'scalar', 'positive'});
            parser.add('rhocp', 4184, 'scalar');
            parser.add('k', [], {'nonnegative', 'scalar'});
            parser.add('Tchwr', 12.2, 'scalar');
            parser.add('Tchws', 5.5, 'scalar');
            parser.add('Tref', [], 'scalar');
            parser.add('Delta', 1, {'scalar', 'positive'});
            parser.add('uvar', 'enthalpy', 'str');
            parser.add('safe', true(), {'scalar', 'bool'});
        end
        args = parser.parse(varargin{:});
        if isempty(args.Tref)
            args.Tref = args.Tchwr;
        end
        if isempty(args.k)
            args.k = 0.005*args.Vmax;
        end
        self = self.update(args, {'Vmax', 'rhocp', 'k', 'Tchwr', 'Tchws', ...
                                  'Tref', 'Delta', 'uvar', 'safe'});
        self = self.init();
    end%function
    
    function [dxdt, u] = model(self, x, u)
        % Uses volume or enthalpy model depending on self.uvar.
        switch self.uvar
        case 'enthalpy'
            [dxdt, u] = self.model_enthalpy(x, u);
        case 'volume'
            [dxdt, u] = self.model_volume(x, u);
        otherwise
            error('Invalid choice for uvar!');
        end
    end%function
    
    function [dxdt, u] = model_volume(self, x, u)
        % Nonlinear dynamic model for storage tank. In terms of volumes and
        % enthalpies. For these purposes, H = rho*C*H*(T - T0) with T
        % temperature. u gives volumetric flow rates (enthalpies are calculated
        % using current parameters).
        %
        % This model is more natural, but the optimizaiton problem models the
        % storage tank using enthalpies.
        narginchk(3, 3);
        
        % Components of x.
        Hcold = x(self.ix.Hcold);
        Vcold = x(self.ix.Vcold);
        Hhot = x(self.ix.Hhot);
        Vhot = x(self.ix.Vhot);
        
        % Components of u.
        vplus = u(self.iu.add);
        hplus = self.rhocp*(self.Tchws - self.Tref);
        vminus = u(self.iu.withdraw);
        hminus = self.rhocp*(self.Tchwr - self.Tref);
        
        % Add safety factor to flows.
        if self.safe
            lamplus = self.safetylam*(1 - Vcold/self.Vmax - self.safetyfrac);
            vplus = vplus*(1 - exp(lamplus));
            
            lamminus = self.safetylam*(Vcold/self.Vmax - self.safetyfrac);
            vminus = vminus*(1 - exp(lamminus));
        end
        
        % Nonlinear model. Need to be careful about division.
        hcold = self.safediv(Hcold, Vcold, hplus);
        hhot = self.safediv(Hhot, Vhot, hminus);
        dxdt = [
            hplus*vplus - hcold*vminus + self.k*(hhot - hcold);
            vplus - vminus;
            -hhot*vplus + hminus*vminus - self.k*(hhot - hcold);
            -vplus + vminus;
        ];
        
        % Also return u, including safety factor.
        u = [vplus; vminus];
    end%function

    function [dxdt, u] = model_enthalpy(self, x, u)
        % Nonlinear dynamic model for storage tank. In terms of volumes and
        % enthalpies. For these purposes, H = rho*C*H*(T - T0) with T
        % temperature. u gives enthalpy added to and withdrawn from the cold
        % section (volumetric flows are calculated using pars).
        narginchk(3, 3);
        hcold_ref = self.rhocp*(self.Tchws - self.Tref);
        
        % Convert u from enthalpy to volume.
        hplus = u(self.iu.add);
        vplus = hplus/hcold_ref;
        hminus = u(self.iu.withdraw);
        Hcold = x(self.ix.Hcold);
        Vcold = x(self.ix.Vcold);
        vminus = hminus*self.safediv(Vcold, Hcold, 1/hcold_ref);
        u = [vplus; vminus];
        
        % Call volume model.
        [dxdt, uvol] = self.model_volume(x, u);
        
        % Convert u from volume back to enthalpy.
        u(self.iu.add) = uvol(self.iu.add)*hcold_ref;
        u(self.iu.withdraw) = uvol(self.iu.withdraw) ...
                              *self.safediv(Hcold, Vcold, hcold_ref);
    end%function
    
    function [x, uget] = simulate(self, x0, u)
        % [x, u] = self.simulate(x0, u)
        %
        % Simulates the storage tank starting from x0 subject to inputs u.
        %
        % Note that all x (including x0) are returned as a matrix. Also, the
        % returned value of u accounts for any reductions due to tank
        % capacity constraints.
        %
        % u can have size [2, Nt] or [1, Nt].
        narginchk(3, 3);
        Nt = size(u, 2);
        x = NaN(self.Nx, Nt + 1);
        uget = NaN(size(u));
        x(:,1) = x0;
        for t = 1:Nt
            try
                [x(:,t + 1), uget(:,t)] = self.step(x(:,t), u(:,t));
            catch err
                warning('Integration failed at step %d: %s', t, err.message);
                break
            end
        end
    end%function
    
    function [x, u, y] = step(self, x, u)
        % [x, u, y] = self.step(x, u)
        %
        % Runs one step of the integrator from states x for inputs u.
        %
        % u can be given as a 2-element vector [add; withdraw], or as a scalar
        % (u(add) - u(withdraw)). u should be in enthalpy or volume depending
        % on self.uvar.
        %
        % The returned value of u gives the actual value of u experienced by
        % the system. If self.safe is true, it may be less than the requested
        % u to avoid violating constraints.
        %
        % The output y is the relevant scalar quantity (cold enthalpy or volume)
        % adjusted to subtract out any safety capacity.
        narginchk(3, 3);
        
        % Check scalar u.
        uscalar = isscalar(u);
        if uscalar
            switch self.uvar
            case 'enthalpy'
                u = [min(u, 0); -max(u, 0)];
            case 'volume'
                u = [max(u, 0); -min(u, 0)];
            end
        end
        
        % Use integrator.
        [x, u] = self.integrator(x, u);
        x = full(x);
        u = full(u);
        
        % Redo scalar u.
        if uscalar
            u = u(1) - u(2);
        end
        
        % Choose y.
        switch self.uvar
        case 'enthalpy'
            i = self.ix.Hcold;
        case 'volume'
            i = self.ix.Vcold;
        end
        y = x(i)*(1 - self.safe*self.safetyfrac);
    end%function
    
    function plot(self, x, lspecs, fig)
        % self.plot(x, [lspecs={'-ob', '-sr'}], [fig])
        %
        % Plots the given tank trajectory. lspecs should give lspecs for
        % plotting the cold and hot states respectively. fig should be the
        % figure handle to use for the plots.
        narginchk(2, 4);
        if nargin() < 3
            lspecs = {'-ob', '-sr'};
        end
        if nargin() < 4
            fig = figure();
        end
        figure(fig);
        t = (0:(size(x, 2) - 1))*self.Delta;
        
        % Plot volume.
        subplot(3, 1, 1);
        hold('on');
        plot(t, x(self.ix.Vcold,:), lspecs{1}, ...
             t, x(self.ix.Vhot,:), lspecs{2});
        ylabel('Volume (m^3)');
        
        % Plot temperature.
        Tcold = self.Tref + x(self.ix.Hcold,:)./(self.rhocp*x(self.ix.Vcold,:));
        Thot = self.Tref + x(self.ix.Hhot,:)./(self.rhocp*x(self.ix.Vhot,:));
        subplot(3, 1, 2);
        hold('on');
        plot(t, Tcold, lspecs{1}, t, Thot, lspecs{2});
        ylabel('Temperature (*C)');
        
        % Plot enthalpy.
        subplot(3, 1, 3);
        hold('on');
        plot(t, x(self.ix.Hcold,:), lspecs{1}, ...
            t, x(self.ix.Hhot,:), lspecs{2});
        ylabel('Enthalpy (kJ)');
        xlabel('Time (h)');
    end%function
end%methods

methods (Access=protected)
    function self = update(self, pars, parnames)
        % Updates self using the fields of pars. parnames defaults to
        % fieldnames(pars).
        narginchk(2, 3);
        if nargin() < 3
            parnames = fieldnames(pars);
        end
        for i = 1:length(parnames)
            p = parnames{i};
            self.(p) = pars.(p);
        end
    end%methods
    
    function self = init(self)
        % Initialize the integrator object.
        narginchk(1, 1);
        
        % Get symbolic expressions and make raw Casadi Integrator object.
        x = casadi.SX.sym('x', self.Nx);
        u = casadi.SX.sym('u', self.Nu);
        [dxdt, uget] = self.model(x, u);
        problem = struct('x', x, 'p', u, 'ode', dxdt, 'quad', uget/self.Delta);
        options = struct('t0', 0, 'tf', self.Delta);
        raw_integrator = casadi.integrator('raw_integrator', 'idas', ...
                                           problem, options);
        
        % Wrap in a nice Casadi.Function.
        x = casadi.MX.sym('x', self.Nx);
        u = casadi.MX.sym('u', self.Nu);
        result = raw_integrator('x0', x, 'p', u);
        self.integrator = casadi.Function('integrator', {x, u}, ...
                                          {result.xf, result.qf});
    end%methods
end%methods

methods (Static)
    function z = safediv(x, y, lim)
        % function z = safediv(x, y, [lim=1])
        %
        % "Safe" division function x/y for positive numbers x and y. lim
        % chooses limit for 0/0.
        narginchk(2, 3);
        if nargin() < 3
            lim = 1;
        end
        z = (x + lim*eps())/(y + eps());
    end%function
end%methods

end%classdef
