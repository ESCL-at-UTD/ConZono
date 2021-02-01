function out = halfspaceIntersection(obj,H)
m = size(H.H,2)-1;
if obj.n ~= m
    disp(['Cannot add a zonotope in ',num2str(obj.n),...
        ' dimensions to a halfspace in ',num2str(m),' dimensions!'])
else
    out = copy(obj);
    for j = 1:size(H.H,1) % considers each of the hyperplane independently
        if abs(H.H(j,end)-H.H(j,1:end-1)*out.c) < sum(abs(H.H(j,1:end-1)*out.G)) % Checks for zonotope intersection with a halfspace
            d_max = H.H(j,end)-H.H(j,1:end-1)*out.c+sum(abs(H.H(j,1:end-1)*out.G)); % Computes the normalized distance of the hyperplane from the farthest vertex of the zonotope
            out.A = [out.A zeros(out.nC,1); H.H(j,1:end-1)*out.G d_max/2]; % Adds constraints
            out.b = [out.b; H.H(j,end)-H.H(j,1:end-1)*out.c-d_max/2];
            out.c = out.c;
            out.G = [out.G, zeros(out.n,1)];
        end
    end
end
end