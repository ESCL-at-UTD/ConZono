function out = volumeratio(obj1,obj2)
if obj1.n ~= obj2.n
	disp(['Cannot compute the Pontryagin difference between a zonotope in ',num2str(obj1.n),...
        ' dimensions and a zonotope in ',num2str(obj2.n),' dimensions!'])
else
	out = (obj1.volume/obj2.volume)^(1/obj1.n);
end