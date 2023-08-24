# -*- coding: utf-8 -*-
module Combustion

#         function [wdot] = wdot_mak(self, z, T, Y, L)
#             % Mass action kinetics methane combustion rate [kg/(m³.s)].
#             k0 = 1.1e+07;
#             Ea = 83680.0;

#             X = self.mass_to_mole_fraction(Y);
#             C = (X * self.PRESSURE ./ (self.GAS_CONSTANT .* T));
#             k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#             rt = k .* C(:, 1) .* C(:, 2).^0.5;

#             wdot = rt * (self.mw .* self.SPECIES_COEFS);
#         endfunction

#         function [wdot] = wdot_ebu(self, z, T, Y, L)
#             % Eddy break-up kinetics methane combustion rate [kg/(m³.s)].
#             cr = 4.000e+00;
#             bo = 4.375e+00;
#             k0 = 1.600e+10;
#             Ea = 1.081e+05;

#             k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#             rho = self.density_mass(T, self.PRESSURE, Y);
#             yf = Y(:, 1);
#             yo = Y(:, 2);

#             % TODO implement this in ProjectData
#             ke = z ./ L;

#             R_ebu = (rho.^1) .* cr .* ke .* min(yf, yo ./ bo);
#             R_arr = (rho.^2) .* yf .* yo .* k;

#             rt = min(R_ebu, R_arr) / self.mw(1);

#             wdot = rt * (self.mw .* self.SPECIES_COEFS);
#         endfunction

#         function [mu] = viscosity(self, T, Y)
#             % Gas molecular viscosity [Pa.s].
#             mu = 1.0e-05 * (0.1672 * sqrt(T) - 1.058);
#         endfunction

#         function [k] = thermal_conductivity(self, T, Y)
#             % Gas thermal conductivity [W/(m³.K)].
#             k = 1.581e-17;
#             k = T .* k - 9.463e-14;
#             k = T .* k + 2.202e-10;
#             k = T .* k - 2.377e-07;
#             k = T .* k + 1.709e-04;
#             k = T .* k - 7.494e-03;
#         endfunction

#     % Species stoichiometric coefficients.
#     SPECIES_COEFS = [-1.0, -2.0, 1.0, 2.0, 0.0, 0.0];

#     % Operating pressure [Pa].
#     PRESSURE = 101325.0;

end # (module Combustion)