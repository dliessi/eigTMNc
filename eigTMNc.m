function [mu, saved] = eigTMNc(model, par, s, h, M, N, varargin)
% eigTMNc computes eigenvalues of the evolution operator T(s + h, s)
%   of linear coupled Renewal Equations and Retarded Functional
%   Differential Equations (RFE/RFDE).
%
%   mu = eigTMNc(model, par, s, h, M, N)
%
%   model: function handle that defines the equations (see the
%          "model_*.m" files for example model definitions);
%          will be called as model(par)
%   par: model parameters
%   M, N: discretization indices
%
%   eigTMNc can reuse some of the computed data from a previous run:
%
%   [mu, saved] = eigTMNc(model, par, s, h, M, N, saved)
%
%   where "saved" is a struct containing the computed data that is
%   independent of the model parameters.
%
%   eigTMNc will ignore the input "saved" variable if it is not a
%   well-formed struct or if it does not correspond to the method
%   parameters for the current run.
%   This can be useful in loops (where reusing data is crucial)
%   since it avoids the need for a special case for the first
%   iteration.
%   The following code will output a warning and recompute the data
%   from scratch in the first iteration, while it will reuse data in
%   the following iterations.
%
%   saved = [];
%   for i = 1 : maxiter
%       [mu, saved] = eigTMNc(model, par, s, h, M, N, saved);
%       % [...]
%   end
%
%   Authors: Davide Liessi <surname.name@spes.uniud.it>
%                and Dimitri Breda
%            Department of Mathematics, Computer Science and Physics
%            University of Udine, Italy
%
%   v0.3 2017-10-23


%% Check sanity of inputs.

assert(isa(model, 'function_handle'), ...
       'model must be a function handle.');
% Sanity of par can be verified by model.
assert(isscalar(s) && isreal(s), 's must be a real scalar.');
assert(isscalar(h) && isreal(h) && h > 0, ...
       'h must be a positive real scalar.');
assert(isscalar(M) && isreal(M) && ceil(M) == floor(M) && M > 0, ...
       'M must be a positive integer.');
assert(isscalar(N) && isreal(N) && ceil(N) == floor(N) && N > 0, ...
       'N must be a positive integer.');
% Postpone sanity check for saved data until we have dX, dY.


%% Model-related definitions.

[dX, dY, tau, AXY, AYY, BXY, BYY, CXX, CXY, CYX, CYY] = model(par);
% tau includes 0 and is sorted.

% Now that we have dX, dY we can check sanity of saved data.
saved_sane = false;
if ~isempty(varargin)
    saved = varargin{1};
    if isstruct(saved)
        saved_sane = ...
            isfield(saved, 'M') && (saved.M == M) && ...
            isfield(saved, 'N') && (saved.N == N) && ...
            isfield(saved, 'dX') && (saved.dX == dX) && ...
            isfield(saved, 'dY') && (saved.dY == dY) && ...
            isfield(saved, 'h') && (saved.h == h);
        if ~saved_sane
            warning('eigTMNc:savedDataMismatch', ...
                    ['Saved data from previous runs ' ...
                     'does not match the other inputs.']);
        end
    else
        warning('eigTMNc:savedDataNotStruct', ...
                ['Saved data from previous runs ' ...
                 'must be passed as a struct.']);
    end
end

TAU = tau(end);
p = length(tau) - 1;

flags.AXY = ~isempty(AXY);
flags.AYY = ~isempty(AYY);
flags.BXY = cell(p, 1);
flags.BYY = cell(p, 1);
flags.CXX = cell(p, 1);
flags.CXY = cell(p, 1);
flags.CYX = cell(p, 1);
flags.CYY = cell(p, 1);
for k = 1 : p
    flags.BXY{k} = ~isempty(BXY{k});
    flags.BYY{k} = ~isempty(BYY{k});
    flags.CXX{k} = ~isempty(CXX{k});
    flags.CXY{k} = ~isempty(CXY{k});
    flags.CYX{k} = ~isempty(CYX{k});
    flags.CYY{k} = ~isempty(CYY{k});
end


%% Method-related definitions.

% Chebyshev zeros in [-1, 1].
chebM = cos((0 : M) * pi / M);
% Chebyshev extrema in [-1, 1].
chebN = cos((2 * (1 : N) - 1) * pi / (2 * N));

% Chebyshev zeros barycentric interpolation weights in [-1, 1].
wbM = [1 / 2, ones(1, M - 1), 1 / 2] .* (- 1) .^ (0 : M);
% Chebyshev extrema barycentric interpolation weights in [-1, 1].
wbN = sin((2 * (1 : N) - 1) * pi / (2 * N)) .* ...
      (- 1) .^ ((1 : N) - 1);

% Clenshaw-Curtis quadrature weights.
wqM = clencurt(M);
wqM_XX = repmat(reshape(wqM, [1 1 M + 1]), [dX dX]);
wqM_XY = repmat(reshape(wqM, [1 1 M + 1]), [dX dY]);
wqM_YX = repmat(reshape(wqM, [1 1 M + 1]), [dY dX]);
wqM_YY = repmat(reshape(wqM, [1 1 M + 1]), [dY dY]);

% Mesh back in time.
if h >= TAU
    OmegaM = (chebM - 1) * TAU / 2;
else
    Q = ceil(TAU / h);
    OmegaM = zeros(Q, M + 1);
    for q = 1 : Q - 1
        OmegaM(q, :) = (chebM - 1) * h / 2 - (q - 1) * h;
    end
    OmegaM(Q, :) = (chebM - 1) * (TAU - (Q - 1) * h) / 2 - ...
                   (Q - 1) * h;
end
% Mesh forward in time.
OmegaN = (1 - chebN) * h / 2;

% Allocate useful variables.
IM = eye(M + 1);
IN = eye(N);
IdX = eye(dX);
IdY = eye(dY);
vq_XX = zeros([dX, dX, M + 1]);
vq_XY = zeros([dX, dY, M + 1]);
vq_YX = zeros([dY, dX, M + 1]);
vq_YY = zeros([dY, dY, M + 1]);


%% The matrix T_{M}^{(1)}.

if saved_sane
    TM1 = saved.TM1;
else
    if h >= TAU
        TM1 = zeros((dX + dY) * (M + 1));
        TM1(dX * (M + 1) + (1 : dY * (M + 1)), ...
            dX * (M + 1) + (1 : dY)) = kron(ones(M + 1, 1), IdY);
        if h == TAU
            TM1(dX * M + (1 : dX), 1 : dX) = IdX;
        end
    else
        TM1 = zeros((dX + dY) * (Q * M + 1));
        TM1(dX * M + 1 : dX * (Q - 1) * M, ...
            1 : dX * (Q - 2) * M) = kron(eye((Q - 2) * M), IdX);
        TM1(dX * (Q * M + 1) + (1 : dY * M), ...
            dX * (Q * M + 1) + (1 : dY)) = kron(ones(M, 1), IdY);
        TM1(dX * (Q * M + 1) + (dY * M + 1 : dY * (Q - 1) * M), ...
            dX * (Q * M + 1) + (1 : dY * (Q - 2) * M)) = ...
            kron(eye((Q - 2) * M), IdY);
        theta = h + OmegaM(Q, :)';
        for m = 0 : M
            lMmq = mybarint(OmegaM(Q - 1, :), wbM, IM(m + 1, :), ...
                            theta);
            TM1(dX * (Q - 1) * M + 1 : dX * (Q * M + 1), ...
                dX * ((Q - 2) * M + m) + (1 : dX)) = ...
                kron(lMmq, IdX);
            TM1(dX * (Q * M + 1) + ...
                (dY * (Q - 1) * M + 1 : dY * (Q * M + 1)), ...
                dX * (Q * M + 1) + ...
                dY * ((Q - 2) * M + m) + (1 : dY)) = ...
                kron(lMmq, IdY);
        end
    end
end


%% The matrix T_{M, N}^{(2)}.

if saved_sane
    TMN2 = saved.TMN2;
else
    if h >= TAU
        TMN2 = zeros((dX + dY) * (M + 1), (dX + dY) * N);
        theta = h + OmegaM';
        for n = 1 : N
            lNn = mybarint(OmegaN, wbN, IN(n, :), theta);
            TMN2(1 : dX * M, dX * (n - 1) + (1 : dX)) = ...
                kron(lNn(1 : M), IdX);
            for m = 0 : M - 1
                int_lNn = ...
                    theta(m + 1) / 2 * ...
                    sum(wqM .* mybarint(OmegaN, wbN, IN(n, :), ...
                        (chebM + 1) * theta(m + 1) / 2));
                TMN2(dX * (M + 1) + dY * m + (1 : dY), ...
                     dX * N + dY * (n - 1) + (1 : dY)) = ...
                    int_lNn * IdY;
            end
            if h > TAU
                TMN2(dX * M + (1 : dX), ...
                     dX * (n - 1) + (1 : dX)) = ...
                    kron(lNn(M + 1), IdX);
                int_lNn = ...
                    theta(M + 1) / 2 * ...
                    sum(wqM .* mybarint(OmegaN, wbN, IN(n, :), ...
                        (chebM + 1) * theta(M + 1) / 2));
                TMN2(dX * (M + 1) + dY * M + (1 : dY), ...
                     dX * N + dY * (n - 1) + (1 : dY)) = ...
                    int_lNn * IdY;
            end
        end
    else
        TMN2 = zeros((dX + dY) * (Q * M + 1), (dX + dY) * N);
        theta = h + OmegaM(1, 1 : M)';
        for n = 1 : N
            lNn = mybarint(OmegaN, wbN, IN(n, :), theta);
            TMN2(1 : dX * M, dX * (n - 1) + (1 : dX)) = ...
                kron(lNn, IdX);
            for m = 0 : M - 1
                int_lNn = ...
                    theta(m + 1) / 2 * ...
                    sum(wqM .* mybarint(OmegaN, wbN, IN(n, :), ...
                        (chebM + 1) * theta(m + 1) / 2));
                TMN2(dX * (Q * M + 1) + dY * m + (1 : dY), ...
                     dX * N + dY * (n - 1) + (1 : dY)) = ...
                    int_lNn * IdY;
            end
        end
    end
end


%% The matrix U_{M, N}^{(1)}.

if h >= TAU
    Nhat = find(OmegaN <= TAU, 1, 'last');
    if isempty(Nhat)
        Nhat = 0;
    end
    GammaXY = cell(N, 1);
    GammaYY = cell(N, 1);
    GammaXY(:) = {zeros(dX, dY)};
    GammaYY(:) = {zeros(dY, dY)};
    for n = 1 : Nhat
        if flags.AXY
            GammaXY{n} = GammaXY{n} + AXY(s + OmegaN(n), par);
        end
        if flags.AYY
            GammaYY{n} = GammaYY{n} + AYY(s + OmegaN(n), par);
        end
        K = find(tau < OmegaN(n), 1, 'last') - 1;
        for k = 1 : K
            if flags.BXY{k}
                GammaXY{n} = GammaXY{n} + ...
                             BXY{k}(s + OmegaN(n), par);
            end
            if flags.BYY{k}
                GammaYY{n} = GammaYY{n} + ...
                             BYY{k}(s + OmegaN(n), par);
            end
            if flags.CXY{k} || flags.CYY{k}
                nqM = (chebM * (tau(k + 1) - tau(k)) - ...
                       tau(k + 1) - tau(k)) / 2;
                int_scale = (tau(k + 1) - tau(k)) / 2;
                vq_XY = zeros([dX, dY, M + 1]);
                vq_YY = zeros([dY, dY, M + 1]);
                if flags.CXY{k}
                    for ii = 1 : M + 1
                        vq_XY(:, :, ii) = ...
                            CXY{k}(s + OmegaN(n), nqM(ii), par);
                    end
                    int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                    GammaXY{n} = GammaXY{n} + int_XY;
                end
                if flags.CYY{k}
                    for ii = 1 : M + 1
                        vq_YY(:, :, ii) = ...
                            CYY{k}(s + OmegaN(n), nqM(ii), par);
                    end
                    int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                    GammaYY{n} = GammaYY{n} + int_YY;
                end
            end
        end
        if flags.CXY{K + 1} || flags.CYY{K + 1}
            nqM = (chebM * (OmegaN(n) - tau(K + 1)) - ...
                   OmegaN(n) - tau(K + 1)) / 2;
            int_scale = (OmegaN(n) - tau(K + 1)) / 2;
            vq_XY = zeros([dX, dY, M + 1]);
            vq_YY = zeros([dY, dY, M + 1]);
            if flags.CXY{K + 1}
                for ii = 1 : M + 1
                    vq_XY(:, :, ii) = ...
                        CXY{K + 1}(s + OmegaN(n), nqM(ii), par);
                end
                int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                GammaXY{n} = GammaXY{n} + int_XY;
            end
            if flags.CYY{K + 1}
                for ii = 1 : M + 1
                    vq_YY(:, :, ii) = ...
                        CYY{K + 1}(s + OmegaN(n), nqM(ii), par);
                end
                int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                GammaYY{n} = GammaYY{n} + int_YY;
            end
        end
    end
    for n = Nhat + 1 : N
        if flags.AXY
            GammaXY{n} = GammaXY{n} + AXY(s + OmegaN(n), par);
        end
        if flags.AYY
            GammaYY{n} = GammaYY{n} + AYY(s + OmegaN(n), par);
        end
        for k = 1 : p
            if flags.BXY{k}
                GammaXY{n} = GammaXY{n} + ...
                             BXY{k}(s + OmegaN(n), par);
            end
            if flags.BYY{k}
                GammaYY{n} = GammaYY{n} + ...
                             BYY{k}(s + OmegaN(n), par);
            end
            if flags.CXY{k} || flags.CYY{k}
                nqM = (chebM * (tau(k + 1) - tau(k)) - ...
                       tau(k + 1) - tau(k)) / 2;
                int_scale = (tau(k + 1) - tau(k)) / 2;
                vq_XY = zeros([dX, dY, M + 1]);
                vq_YY = zeros([dY, dY, M + 1]);
                if flags.CXY{k}
                    for ii = 1 : M + 1
                        vq_XY(:, :, ii) = ...
                            CXY{k}(s + OmegaN(n), nqM(ii), par);
                    end
                    int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                    GammaXY{n} = GammaXY{n} + int_XY;
                end
                if flags.CYY{k}
                    for ii = 1 : M + 1
                        vq_YY(:, :, ii) = ...
                            CYY{k}(s + OmegaN(n), nqM(ii), par);
                    end
                    int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                    GammaYY{n} = GammaYY{n} + int_YY;
                end
            end
        end
    end
    ThetaXX = cell(N, M + 1);
    ThetaXY = cell(N, M + 1);
    ThetaYX = cell(N, M + 1);
    ThetaYY = cell(N, M + 1);
    ThetaXX(:, :) = {zeros(dX, dX)};
    ThetaXY(:, :) = {zeros(dX, dY)};
    ThetaYX(:, :) = {zeros(dY, dX)};
    ThetaYY(:, :) = {zeros(dY, dY)};
    for m = 0 : M
        for n = 1 : Nhat
            K = find(tau < OmegaN(n), 1, 'last') - 1;
            for k = K + 1 : p
                if flags.BXY{k} || flags.BYY{k}
                    lMm = mybarint(OmegaM, wbM, IM(m + 1, :), ...
                                   OmegaN(n) - tau(k + 1));
                    if flags.BXY{k}
                        ThetaXY{n, m + 1} = ...
                            ThetaXY{n, m + 1} + ...
                            BXY{k}(s + OmegaN(n), par) * lMm;
                    end
                    if flags.BYY{k}
                        ThetaYY{n, m + 1} = ...
                            ThetaYY{n, m + 1} + ...
                            BYY{k}(s + OmegaN(n), par) * lMm;
                    end
                end
            end
            if flags.CXX{K + 1} || flags.CXY{K + 1} || ...
               flags.CYX{K + 1} || flags.CYY{K + 1}
                nqM = (chebM * (tau(K + 2) - OmegaN(n)) - ...
                       tau(K + 2) - OmegaN(n)) / 2;
                int_scale = (tau(K + 2) - OmegaN(n)) / 2;
                vq_lMm = mybarint(OmegaM, wbM, IM(m + 1, :), ...
                                  OmegaN(n) + nqM);
                if flags.CXX{K + 1}
                    for ii = 1 : M + 1
                        vq_XX(:, :, ii) = ...
                            CXX{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                    end
                    int_XX = int_scale * sum(wqM_XX .* vq_XX, 3);
                    ThetaXX{n, m + 1} = ThetaXX{n, m + 1} + int_XX;
                end
                if flags.CXY{K + 1}
                    for ii = 1 : M + 1
                        vq_XY(:, :, ii) = ...
                            CXY{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                    end
                    int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                    ThetaXY{n, m + 1} = ThetaXY{n, m + 1} + int_XY;
                end
                if flags.CYX{K + 1}
                    for ii = 1 : M + 1
                        vq_YX(:, :, ii) = ...
                            CYX{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                    end
                    int_YX = int_scale * sum(wqM_YX .* vq_YX, 3);
                    ThetaYX{n, m + 1} = ThetaYX{n, m + 1} + int_YX;
                end
                if flags.CYY{K + 1}
                    for ii = 1 : M + 1
                        vq_YY(:, :, ii) = ...
                            CYY{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                    end
                    int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                    ThetaYY{n, m + 1} = ThetaYY{n, m + 1} + int_YY;
                end
            end
            for k = K + 2 : p
                if flags.CXX{k} || flags.CXY{k} || ...
                   flags.CYX{k} || flags.CYY{k}
                    nqM = (chebM * (tau(k + 1) - tau(k)) - ...
                           tau(k + 1) - tau(k)) / 2;
                    int_scale = (tau(k + 1) - tau(k)) / 2;
                    vq_lMm = mybarint(OmegaM, wbM, IM(m + 1, :), ...
                                      OmegaN(n) + nqM);
                    if flags.CXX{k}
                        for ii = 1 : M + 1
                            vq_XX(:, :, ii) = ...
                                CXX{k}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                        end
                        int_XX = int_scale * ...
                                 sum(wqM_XX .* vq_XX, 3);
                        ThetaXX{n, m + 1} = ThetaXX{n, m + 1} + ...
                                            int_XX;
                    end
                    if flags.CXY{k}
                        for ii = 1 : M + 1
                            vq_XY(:, :, ii) = ...
                                CXY{k}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                        end
                        int_XY = int_scale * ...
                                 sum(wqM_XY .* vq_XY, 3);
                        ThetaXY{n, m + 1} = ThetaXY{n, m + 1} + ...
                                            int_XY;
                    end
                    if flags.CYX{k}
                        for ii = 1 : M + 1
                            vq_YX(:, :, ii) = ...
                                CYX{k}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                        end
                        int_YX = int_scale * ...
                                 sum(wqM_YX .* vq_YX, 3);
                        ThetaYX{n, m + 1} = ThetaYX{n, m + 1} + ...
                                            int_YX;
                    end
                    if flags.CYY{k}
                        for ii = 1 : M + 1
                            vq_YY(:, :, ii) = ...
                                CYY{k}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lMm(ii);
                        end
                        int_YY = int_scale * ...
                                 sum(wqM_YY .* vq_YY, 3);
                        ThetaYY{n, m + 1} = ThetaYY{n, m + 1} + ...
                                            int_YY;
                    end
                end
            end
        end
    end
    UMN1 = zeros((dX + dY) * N, (dX + dY) * (M + 1));
    for m = 0 : M
        for n = 1 : Nhat
            UMN1(dX * (n - 1) + (1 : dX), ...
                 dX * m + (1 : dX)) = ...
                ThetaXX{n, m + 1};
            UMN1(dX * (n - 1) + (1 : dX), ...
                 dX * (M + 1) + dY * m + (1 : dY)) = ...
                ThetaXY{n, m + 1};
            UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                 dX * m + (1 : dX)) = ...
                ThetaYX{n, m + 1};
            UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                 dX * (M + 1) + dY * m + (1 : dY)) = ...
                ThetaYY{n, m + 1};
        end
    end
    for n = 1 : N
        UMN1(dX * (n - 1) + (1 : dX), ...
             dX * (M + 1) + (1 : dY)) = ...
            UMN1(dX * (n - 1) + (1 : dX), ...
                 dX * (M + 1) + (1 : dY)) + GammaXY{n};
        UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
             dX * (M + 1) + (1 : dY)) = ...
            UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                 dX * (M + 1) + (1 : dY)) + GammaYY{n};
    end
else
    tNnq = h * repmat(0 : (Q - 1), N, 1) + repmat(OmegaN', 1, Q);
    GammaXY = cell(N, 1);
    GammaYY = cell(N, 1);
    GammaXY(:) = {zeros(dX, dY)};
    GammaYY(:) = {zeros(dY, dY)};
    for n = 1 : N
        if flags.AXY
            GammaXY{n} = AXY(s + OmegaN(n), par);
        end
        if flags.AYY
            GammaYY{n} = AYY(s + OmegaN(n), par);
        end
        K = find(tau < OmegaN(n), 1, 'last') - 1;
        for k = 1 : K
            if flags.BXY{k}
                GammaXY{n} = GammaXY{n} + ...
                             BXY{k}(s + OmegaN(n), par);
            end
            if flags.BYY{k}
                GammaYY{n} = GammaYY{n} + ...
                             BYY{k}(s + OmegaN(n), par);
            end
            if flags.CXY{k} || flags.CYY{k}
                nqM = (chebM * (tau(k + 1) - tau(k)) - ...
                       tau(k + 1) - tau(k)) / 2;
                int_scale = (tau(k + 1) - tau(k)) / 2;
                vq_XY = zeros([dX, dY, M + 1]);
                vq_YY = zeros([dY, dY, M + 1]);
                if flags.CXY{k}
                    for ii = 1 : M + 1
                        vq_XY(:, :, ii) = CXY{k}(s + OmegaN(n), ...
                                                 nqM(ii), par);
                    end
                    int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                    GammaXY{n} = GammaXY{n} + int_XY;
                end
                if flags.CYY{k}
                    for ii = 1 : M + 1
                        vq_YY(:, :, ii) = CYY{k}(s + OmegaN(n), ...
                                                 nqM(ii), par);
                    end
                    int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                    GammaYY{n} = GammaYY{n} + int_YY;
                end
            end
        end
        if flags.CXY{K + 1} || flags.CYY{K + 1}
            nqM = (chebM * (OmegaN(n) - tau(K + 1)) - ...
                   OmegaN(n) - tau(K + 1)) / 2;
            int_scale = (OmegaN(n) - tau(K + 1)) / 2;
            vq_XY = zeros([dX, dY, M + 1]);
            vq_YY = zeros([dY, dY, M + 1]);
            if flags.CXY{K + 1}
                for ii = 1 : M + 1
                    vq_XY(:, :, ii) = CXY{K + 1}(s + OmegaN(n), ...
                                                 nqM(ii), par);
                end
                int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                GammaXY{n} = GammaXY{n} + int_XY;
            end
            if flags.CYY{K + 1}
                for ii = 1 : M + 1
                    vq_YY(:, :, ii) = CYY{K + 1}(s + OmegaN(n), ...
                                                 nqM(ii), par);
                end
                int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                GammaYY{n} = GammaYY{n} + int_YY;
            end
        end
    end
    ThetaXX = cell(N, M + 1, Q);
    ThetaXY = cell(N, M + 1, Q);
    ThetaYX = cell(N, M + 1, Q);
    ThetaYY = cell(N, M + 1, Q);
    ThetaXX(:, :, :) = {zeros(dX, dX)};
    ThetaXY(:, :, :) = {zeros(dX, dY)};
    ThetaYX(:, :, :) = {zeros(dY, dX)};
    ThetaYY(:, :, :) = {zeros(dY, dY)};
    for n = 1 : N
        for q = 1 : Q
            if q ~= Q
                KtNnq = find(tau < tNnq(n, q + 1), 1, 'last') - 1;
                KBnq = min([KtNnq, p]);
                KCnq = min([KtNnq + 1, p]);
            else
                KBnq = p;
                KCnq = p;
            end
            KtNnqm1 = find(tau < tNnq(n, q), 1, 'last') - 1;
            for k = KtNnqm1 + 1 : KBnq
                if flags.BXY{k} || flags.BYY{k}
                    for m = 0 : M
                        lMmq = mybarint(OmegaM(q, :), wbM, ...
                                        IM(m + 1, :), ...
                                        OmegaN(n) - tau(k + 1));
                        if flags.BXY{k}
                            ThetaXY{n, m + 1, q} = ...
                                ThetaXY{n, m + 1, q} + ...
                                BXY{k}(s + OmegaN(n), par) * lMmq;
                        end
                        if flags.BYY{k}
                            ThetaYY{n, m + 1, q} = ...
                                ThetaYY{n, m + 1, q} + ...
                                BYY{k}(s + OmegaN(n), par) * lMmq;
                        end
                    end
                end
            end
            for k = KtNnqm1 + 1 : KCnq
                if flags.CXX{k} || flags.CXY{k} || ...
                   flags.CYX{k} || flags.CYY{k}
                    if q ~= Q
                        akq = max([-tau(k + 1), -tNnq(n, q + 1)]);
                    else
                        akq = -tau(k + 1);
                    end
                    bkq = min([-tau(k), -tNnq(n, q)]);
                    if akq ~= bkq
                        nqM = (chebM * (bkq - akq) + akq + bkq) / 2;
                        int_scale = (bkq - akq) / 2;
                        for m = 0 : M
                            vq_lMmq = ...
                                mybarint(OmegaM(q, :), wbM, ...
                                         IM(m + 1, :), ...
                                         OmegaN(n) + nqM);
                            if flags.CXX{k}
                                for ii = 1 : M + 1
                                    vq_XX(:, :, ii) = ...
                                        CXX{k}(s + OmegaN(n), ...
                                               nqM(ii), par) * ...
                                        vq_lMmq(ii);
                                end
                                int_XX = int_scale * ...
                                         sum(wqM_XX .* vq_XX, 3);
                                ThetaXX{n, m + 1, q} = ...
                                    ThetaXX{n, m + 1, q} + int_XX;
                            end
                            if flags.CXY{k}
                                for ii = 1 : M + 1
                                    vq_XY(:, :, ii) = ...
                                        CXY{k}(s + OmegaN(n), ...
                                               nqM(ii), par) * ...
                                        vq_lMmq(ii);
                                end
                                int_XY = int_scale * ...
                                         sum(wqM_XY .* vq_XY, 3);
                                ThetaXY{n, m + 1, q} = ...
                                    ThetaXY{n, m + 1, q} + int_XY;
                            end
                            if flags.CYX{k}
                                for ii = 1 : M + 1
                                    vq_YX(:, :, ii) = ...
                                        CYX{k}(s + OmegaN(n), ...
                                               nqM(ii), par) * ...
                                        vq_lMmq(ii);
                                end
                                int_YX = int_scale * ...
                                         sum(wqM_YX .* vq_YX, 3);
                                ThetaYX{n, m + 1, q} = ...
                                    ThetaYX{n, m + 1, q} + int_YX;
                            end
                            if flags.CYY{k}
                                for ii = 1 : M + 1
                                    vq_YY(:, :, ii) = ...
                                        CYY{k}(s + OmegaN(n), ...
                                               nqM(ii), par) * ...
                                        vq_lMmq(ii);
                                end
                                int_YY = int_scale * ...
                                         sum(wqM_YY .* vq_YY, 3);
                                ThetaYY{n, m + 1, q} = ...
                                    ThetaYY{n, m + 1, q} + int_YY;
                            end
                        end
                    end
                end
            end
        end
    end
    UMN1 = zeros((dX + dY) * N, (dX + dY) * (Q * M + 1));
    for m = 0 : M
        for q = 1 : Q
            for n = 1 : N
                UMN1(dX * (n - 1) + (1 : dX), ...
                     dX * (m + (q - 1) * M) + (1 : dX)) = ...
                    UMN1(dX * (n - 1) + (1 : dX), ...
                         dX * (m + (q - 1) * M) + (1 : dX)) + ...
                    ThetaXX{n, m + 1, q};
                UMN1(dX * (n - 1) + (1 : dX), ...
                     dX * (Q * M + 1) + ...
                     dY * (m + (q - 1) * M) + (1 : dY)) = ...
                    UMN1(dX * (n - 1) + (1 : dX), ...
                         dX * (Q * M + 1) + ...
                         dY * (m + (q - 1) * M) + (1 : dY)) + ...
                    ThetaXY{n, m + 1, q};
                UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                     dX * (m + (q - 1) * M) + (1 : dX)) = ...
                    UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                         dX * (m + (q - 1) * M) + (1 : dX)) + ...
                    ThetaYX{n, m + 1, q};
                UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                     dX * (Q * M + 1) + ...
                     dY * (m + (q - 1) * M) + (1 : dY)) = ...
                    UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                         dX * (Q * M + 1) + ...
                         dY * (m + (q - 1) * M) + (1 : dY)) + ...
                    ThetaYY{n, m + 1, q};
            end
        end
    end
    for n = 1 : N
        UMN1(dX * (n - 1) + (1 : dX), ...
             dX * (Q * M + 1) + (1 : dY)) = ...
            UMN1(dX * (n - 1) + (1 : dX), ...
                 dX * (Q * M + 1) + (1 : dY)) + ...
            GammaXY{n};
        UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
             dX * (Q * M + 1) + (1 : dY)) = ...
            UMN1(dX * N + dY * (n - 1) + (1 : dY), ...
                 dX * (Q * M + 1) + (1 : dY)) + ...
            GammaYY{n};
    end
end


%% The matrix U_{N}^{(2)}.

% lNn = @(n, t) mybarint(OmegaN, wbN, IN(n, :), t);

LambdaXX = cell(N, N);
LambdaXY = cell(N, N);
LambdaYX = cell(N, N);
LambdaYY = cell(N, N);
LambdaXX(:, :) = {zeros(dX, dX)};
LambdaXY(:, :) = {zeros(dX, dY)};
LambdaYX(:, :) = {zeros(dY, dX)};
LambdaYY(:, :) = {zeros(dY, dY)};
UN2 = zeros((dX + dY) * N);
for n = 1 : N
    for i = 1 : N
        if flags.AXY || flags.AYY
            int_lNn = OmegaN(n) / 2 * ...
                      sum(wqM .* mybarint(OmegaN, wbN, IN(i, :), ...
                          (chebM + 1) * OmegaN(n) / 2));
            if flags.AXY
                LambdaXY{n, i} = LambdaXY{n, i} + ...
                                 AXY(s + OmegaN(n), par) * int_lNn;
            end
            if flags.AYY
                LambdaYY{n, i} = LambdaYY{n, i} + ...
                                 AYY(s + OmegaN(n), par) * int_lNn;
            end
        end
        K = find(tau < OmegaN(n), 1, 'last') - 1;
        for k = 1 : K
            if flags.BXY{k} || flags.BYY{k}
                int_lNn = (OmegaN(n) - tau(k + 1)) / 2 * ...
                          sum(wqM .* ...
                              mybarint(OmegaN, wbN, IN(i, :), ...
                              (chebM + 1) * ...
                              (OmegaN(n) - tau(k + 1)) / 2));
                if flags.BXY{k}
                    LambdaXY{n, i} = ...
                        LambdaXY{n, i} + ...
                        BXY{k}(s + OmegaN(n), par) * int_lNn;
                end
                if flags.BYY{k}
                    LambdaYY{n, i} = ...
                        LambdaYY{n, i} + ...
                        BYY{k}(s + OmegaN(n), par) * int_lNn;
                end
            end
        end
        if K < p
            if flags.CXX{K + 1} || flags.CXY{K + 1} || ...
               flags.CYX{K + 1} || flags.CYY{K + 1}
                nqM = (chebM * (OmegaN(n) - tau(K + 1)) - ...
                       OmegaN(n) - tau(K + 1)) / 2;
                int_scale = (OmegaN(n) - tau(K + 1)) / 2;
                vq_lNn = mybarint(OmegaN, wbN, IN(i, :), ...
                                  OmegaN(n) + nqM);
                vq_int_lNn = zeros(M + 1, 1);
                if flags.CXY{K + 1} || flags.CYY{K + 1}
                    for ii = 1 : M + 1
                        vq_int_lNn(ii) = ...
                            (OmegaN(n) + nqM(ii)) / 2 * ...
                            sum(wqM .* ...
                                mybarint(OmegaN, wbN, IN(i, :), ...
                                         (chebM + 1) * ...
                                         (OmegaN(n) + ...
                                          nqM(ii)) / 2));
                    end
                end
                if flags.CXX{K + 1}
                    for ii = 1 : M + 1
                        vq_XX(:, :, ii) = ...
                            CXX{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lNn(ii);
                    end
                    int_XX = int_scale * sum(wqM_XX .* vq_XX, 3);
                    LambdaXX{n, i} = LambdaXX{n, i} + int_XX;
                end
                if flags.CXY{K + 1}
                    for ii = 1 : M + 1
                        vq_XY(:, :, ii) = ...
                            CXY{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * ...
                            vq_int_lNn(ii);
                    end
                    int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                    LambdaXY{n, i} = LambdaXY{n, i} + int_XY;
                end
                if flags.CYX{K + 1}
                    for ii = 1 : M + 1
                        vq_YX(:, :, ii) = ...
                            CYX{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * vq_lNn(ii);
                    end
                    int_YX = int_scale * sum(wqM_YX .* vq_YX, 3);
                    LambdaYX{n, i} = LambdaYX{n, i} + int_YX;
                end
                if flags.CYY{K + 1}
                    for ii = 1 : M + 1
                        vq_YY(:, :, ii) = ...
                            CYY{K + 1}(s + OmegaN(n), ...
                                       nqM(ii), par) * ...
                            vq_int_lNn(ii);
                    end
                    int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                    LambdaYY{n, i} = LambdaYY{n, i} + int_YY;
                end
            end
        end
        for k = 1 : K
            if flags.CXX{k} || flags.CXY{k} || flags.CYX{k} || ...
               flags.CYY{k}
                nqM = (chebM * (tau(k + 1) - tau(k)) - ...
                       tau(k + 1) - tau(k)) / 2;
                int_scale = (tau(k + 1) - tau(k)) / 2;
                vq_lNn = mybarint(OmegaN, wbN, IN(i, :), ...
                                  OmegaN(n) + nqM);
                vq_int_lNn = zeros(M + 1, 1);
                if flags.CXY{k} || flags.CYY{k}
                    for ii = 1 : M + 1
                        vq_int_lNn(ii) = ...
                            (OmegaN(n) + nqM(ii)) / 2 * ...
                            sum(wqM .* ...
                                mybarint(OmegaN, wbN, IN(i, :), ...
                                         (chebM + 1) * ...
                                         (OmegaN(n) + ...
                                          nqM(ii)) / 2));
                    end
                end
                if flags.CXX{k}
                    for ii = 1 : M + 1
                        vq_XX(:, :, ii) = ...
                            CXX{k}(s + OmegaN(n), ...
                                   nqM(ii), par) * vq_lNn(ii);
                    end
                    int_XX = int_scale * sum(wqM_XX .* vq_XX, 3);
                    LambdaXX{n, i} = LambdaXX{n, i} + int_XX;
                end
                if flags.CXY{k}
                    for ii = 1 : M + 1
                        vq_XY(:, :, ii) = ...
                            CXY{k}(s + OmegaN(n), ...
                                   nqM(ii), par) * vq_int_lNn(ii);
                    end
                    int_XY = int_scale * sum(wqM_XY .* vq_XY, 3);
                    LambdaXY{n, i} = LambdaXY{n, i} + int_XY;
                end
                if flags.CYX{k}
                    for ii = 1 : M + 1
                        vq_YX(:, :, ii) = ...
                            CYX{k}(s + OmegaN(n), ...
                                   nqM(ii), par) * vq_lNn(ii);
                    end
                    int_YX = int_scale * sum(wqM_YX .* vq_YX, 3);
                    LambdaYX{n, i} = LambdaYX{n, i} + int_YX;
                end
                if flags.CYY{k}
                    for ii = 1 : M + 1
                        vq_YY(:, :, ii) = ...
                            CYY{k}(s + OmegaN(n), ...
                                   nqM(ii), par) * vq_int_lNn(ii);
                    end
                    int_YY = int_scale * sum(wqM_YY .* vq_YY, 3);
                    LambdaYY{n, i} = LambdaYY{n, i} + int_YY;
                end
            end
        end
        UN2(dX * (n - 1) + (1 : dX), ...
            dX * (i - 1) + (1 : dX)) = LambdaXX{n, i};
        UN2(dX * (n - 1) + (1 : dX), ...
            dX * N + dY * (i - 1) + (1 : dY)) = LambdaXY{n, i};
        UN2(dX * N + dY * (n - 1) + (1 : dY), ...
            dX * (i - 1) + (1 : dX)) = LambdaYX{n, i};
        UN2(dX * N + dY * (n - 1) + (1 : dY), ...
            dX * N + dY * (i - 1) + (1 : dY)) = LambdaYY{n, i};
    end
end


%% The matrix T_{M, N}.

TMN = TM1 + TMN2 * ((eye((dX + dY) * N) - UN2) \ UMN1);


%% Compute eigenvalues.

mu = eig(TMN);
% Sort mu by decreasing magnitude.
% Work around missing 'ComparisonMethod' option of sort in Octave.
if isreal(mu)
    [mu, ii] = sort([1i; mu], 'descend');
    mu(find(ii == 1)) = [];
else
    mu = sort(mu, 'descend');
end

if ~saved_sane
    saved.M = M;
    saved.N = N;
    saved.dX = dX;
    saved.dY = dY;
    saved.h = h;
    saved.TM1 = TM1;
    saved.TMN2 = TMN2;
end

end


%% Auxiliary functions.

function [w, x] = clencurt(N)
%CLENCURT weights and nodes of Clenshaw-Curtis quadrature
%  [w, x] = CLENCURT(N) returns the weights of the Clenshaw-Curtis
%  quadrature, i.e. the pseudospectral method on N + 1 Chebyshev
%  nodes x in [-1, 1].
%
%  Reference:
%    L. N. Trefethen,
%    Spectral Methods in Matlab,
%    SIAM, 2000,
%    DOI:10.1137/1.9780898719598 .

p = pi * (0 : N)' / N;
x = cos(p);
w = zeros(1, N + 1);
ii = 2 : N;
v = ones(N - 1, 1);
if mod(N, 2) == 0
    w(1) = 1 / (N ^ 2 - 1);
    w(N + 1) = w(1);
    for k = 1 : N / 2 - 1
        v = v - 2 * cos(2 * k * p(ii)) / (4 * k ^ 2 - 1);
    end
    v = v - cos(N * p(ii)) / (N ^ 2 - 1);
else
    w(1) = 1 / N ^ 2;
    w(N + 1) = w(1);
    for	k = 1 : (N - 1) / 2
        v = v - 2 * cos(2 * k * p(ii)) / (4 * k ^ 2 - 1);
    end
end
w(ii) = 2 * v / N;
end

function ff = mybarint(x, b, f, xx)
%MYBARINT Barycentric interpolation
%  Compute the values ff of a function on xx using the barycentric
%  interpolation formula with x interpolation nodes, b barycentric
%  weights and f values of the function on x.
%
%  ff = MYBARINT(x, b, f, xx)
%
%  Reference:
%    J.-P. Berrut and L. N. Trefethen,
%    Barycentric Lagrange Interpolation,
%    SIAM Review, 46(3):501-517, 2004,
%    DOI:10.1137/S0036144502417715 .

n = length(x);

numer = zeros(size(xx));
denom = zeros(size(xx));
exact = zeros(size(xx));
for j = 1 : n
    tdiff = xx - x(j);
    temp = b(j) ./ tdiff;
    numer = numer + temp * f(j);
    denom = denom + temp;
    exact(tdiff == 0) = j;
end
jj = find(exact);
ff = numer ./ denom;
ff(jj) = f(exact(jj));
end
