function [Phi,Psi] = stitch_mrdmd(nodes,level)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% start = min(extractfield(nodes,'start'));
% stop = max(extractfield(nodes,'stop'));
start = min([nodes.start]);
stop = max([nodes.stop]);
t = stop - start + 1;

% levels = extractfield(nodes,'level');
levels = [nodes.level];

nodes = nodes(levels==level);
if(~isempty(nodes))
    % Phi_all = extractfield(nodes,'Phi');
    Phi_all = [nodes.Phi];

    % nmodes = sum(extractfield(nodes,'n'));
    nmodes = sum([nodes.n]);
    Phi = reshape(Phi_all,[size(nodes(1).Phi,1),nmodes]);
    Psi = zeros(nmodes,t);
    
    j=1;
    for i=1:length(nodes)
        n = nodes(i);
        nmodes_i = size(n.Psi,1);
        Psi(j:j+nmodes_i-1,n.start:n.stop) = n.Psi;
        j = j + nmodes_i;
    end
else
    Phi = 0;
    Psi = 0;
end
end