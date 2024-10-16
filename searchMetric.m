% function [sol, prog] = searchMetric(z, vv, W, x, prog_base)
function [sol, prog, cost_vec] = searchMetric(pnum, dz, xop, dHdx, Gamma_fun, Omega_fun, x, Jac, prog_base, proj_flag)
% dz = [[z, H + S]] = {z, H + S} + (z, H + S)
% dz = (Omega + Gamma) * (dHdx + dSdx)
cost_vec = [];

parfor k = 1:size(xop, 1)
    Gk = subs(Gamma_fun, x, xop(k, :)'); % As shown in Gibbs free eenrgy, we want the free energy to decrease
    Ok = subs(Omega_fun, x, xop(k, :)');

    Jk = subs(Jac, x, xop(k, :)');
    dSdx = Jk' * pnum;

    % consider metric projection
    if proj_flag
        
        rhs = Ok * (dHdx(k, :)' + dSdx);
        dSdx_proj = dSdx - dSdx' * dHdx(k, :)' / (dHdx(k, :) * dHdx(k, :)') * dHdx(k, :)';
        rhs = rhs + Gk * dSdx_proj;

        cost_vec = [cost_vec; dz(k, :)' - rhs ];
    else
        % do not consider metric projection
        cost_vec = [cost_vec; dz(k, :)' - (Ok + Gk) * (dHdx(k, :)' + dSdx) ];
    end
end
[prog, t] = prog_base.newPos(length(cost_vec)); % 1-norm
% [prog, t] = prog_base.newPos(1); % inf-norm
prog = prog.withPos(t - cost_vec);
prog = prog.withPos(t + cost_vec);

options = spot_sdp_default_options();
options.verbose = 12;
[sol, prog] = prog.minimize(sum(t) * 1e4, @spot_mosek, options);
end

%%
% function [sol, prog] = searchMetric(z, vv, W, x, prog_base)
% % minimize the one-norm |z - A(x)vv|
% cost_vec = [];
% tic
% parfor k = 1:size(vv, 1)
%     Wk = subs(W, x, vv(k, :)');
%     cost_vec = [cost_vec; z(k, :)' - Wk * vv(k, :)'];
% end
% toc
% [prog, t] = prog_base.newPos(length(cost_vec));
% prog = prog.withPos(t - cost_vec);
% prog = prog.withPos(t + cost_vec);
%
% options = spot_sdp_default_options();
% options.verbose = 12;
% [sol, prog] = prog.minimize(sum(t) / 1e2, @spot_mosek, options);
% end