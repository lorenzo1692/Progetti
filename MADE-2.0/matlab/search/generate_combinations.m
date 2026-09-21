function combT = generate_combinations(env)
%GENERATE_COMBINATIONS Equal-turns-per-grade turns/layers layouts.
%
%   combT = GENERATE_COMBINATIONS(env) returns a cell array, one entry per
%   candidate number of layers in env.layers_comb; each entry is a matrix
%   whose rows are candidate turns-per-grade layouts (same turns count
%   repeated across the grades of that layout) whose total turn count
%   falls within env.n_spire_fsbl.

combT = cell(1, numel(env.layers_comb));
for i = 1:numel(env.layers_comb)
    n_layers_i = env.layers_comb(i);
    t = env.turns_comb' .* ones(numel(env.turns_comb), n_layers_i);
    valid = sum(t, 2) > min(env.n_spire_fsbl) & sum(t, 2) < max(env.n_spire_fsbl);
    x = valid .* t;
    x(sum(x, 2) == 0, :) = [];
    combT{i} = x;
end
end
