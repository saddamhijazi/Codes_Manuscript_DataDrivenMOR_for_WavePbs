function nodes = mrdmd(D,level,bin_num,offset,max_levels,max_cycles,do_svht)
% Compute the multi-resolution DMD on the dataset `D`
% 4 times nyquist limit to capture cycles
nyq = 8 * max_cycles;

nodes = [];

% time bin size
bin_size = size(D,2);
if bin_size < nyq
    return

end

% extract subsamples 
step = floor(bin_size/nyq); % max step size to capture cycles
D_sub = D(:,1:step:end);
X = D_sub(:,1:end-1);
Y = D_sub(:,2:end);

% determine rank-reduction
if do_svht
    sv_sub = svd(D_sub);
    tau = svht(D_sub, sv_sub);
    r = sum(sv_sub > tau);
else
    r = min(size(X));
end

% compute dmd
[lambda,Phi] = dmd(X, Y, r);

% frequency cutoff (oscillations per timestep)
rho = max_cycles/bin_size;

% consolidate slow eigenvalues (as boolean mask)
slow = find((abs(log(lambda) / (2 * pi * step))) <= rho);
n = length(slow); % number of slow modes

% extract slow modes (perhaps empty)
lambda = lambda(slow);
Phi = Phi(:,slow);

if(n > 0)

    % vars for the objective function for D (before subsampling)
    Vand = vanderm(power(lambda,1/step),bin_size);
    P = (Phi'*Phi).*conj(Vand*Vand');
    q = conj(diag((Vand*D')*Phi));
    
    % find optimal b solution
    Pl = chol(P,'lower');
    b_opt = (Pl')\(Pl\q); % Optimal vector of amplitudes b
    
    % time evolution
    Psi = diag(b_opt)*Vand;

else

    % zero time evolution
    b_opt = [];
    Psi = zeros([0, bin_size]);

end

% dmd reconstruction
D_dmd = Phi*Psi; 

% remove influence of slow modes
D = D - D_dmd;

% record keeping
% node = type('Node', (object,), {})()
node.level = level;            % level of recursion
node.bin_num = bin_num;        % time bin number
node.bin_size = bin_size;      % time bin size
node.start = offset + 1;           % starting index
node.stop = offset + bin_size; % stopping index
node.step = step;              % step size
node.rho = rho;                % frequency cutoff
node.r = r;                    % rank-reduction
node.n = n;                    % number of extracted modes
node.lambda = lambda;          % extracted eigenvalues
node.Phi = Phi;                % extracted DMD modes
node.Psi = Psi;                % extracted time evolution
node.b_opt = b_opt;            % extracted optimal b vector
nodes = node;

% split data into two and do recursion
if level < max_levels
    split = ceil(bin_size / 2); % where to split
    nodes = [nodes mrdmd(...
        D(:,1:split),...
        level+1,...
        2*bin_num,...
        offset,...
        max_levels,...
        max_cycles,...
        do_svht...
        )];
    nodes = [nodes  mrdmd(...
        D(:,split+1:end),...
        level+1,...
        2*bin_num+1,...
        offset+split,...
        max_levels,...
        max_cycles,...
        do_svht...
        )];
end