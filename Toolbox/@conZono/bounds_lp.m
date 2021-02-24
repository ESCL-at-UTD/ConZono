function out = bounds_lp(obj)

ng = obj.nG;

for i = 1:ng
    s = zeros(ng,1);
    s(i) = 1;
    lp = Opt('f',s','Ae',obj.A,'be',obj.b,'lb',-ones(ng,1),'ub',ones(ng,1)); % Formulates the linear program
    opt = mpt_solve(lp); % Solves the LP.
    out(i,1) = s'*opt.xopt;
end

for i = 1:ng
    s = zeros(ng,1);
    s(i) = 1;
    lp = Opt('f',-s','Ae',obj.A,'be',obj.b,'lb',-ones(ng,1),'ub',ones(ng,1)); % Formulates the linear program
    opt = mpt_solve(lp); % Solves the LP.
    out(i,2) = s'*opt.xopt;
end

end