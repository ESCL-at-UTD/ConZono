function out = plus(obj1,obj2)
arguments
    obj1     Zono
    obj2     
end
checkObj2(obj1,obj2);

if isa(obj2,'Zono')
    out = conZono;
    out.c = obj1.c + obj2.c;
    out.G = [obj1.G obj2.G];
else %obj2 is a column vector
    out = copy(obj1);
    out.c = out.c + obj2; %we simply translate the original conZono
end
end


function checkObj2(obj1, obj2)
if isa(obj2,'Zono')
    if obj1.n ~= obj2.n
        eid = 'conZono:SizesNotEqual';
        msg = ['Cannot add a zonotope in ',num2str(obj1.n),...
            ' dimensions to a zonotope in ',num2str(obj2.n),' dimensions!'];
        throwAsCaller(MException(eid,msg))
    end
elseif isnumeric(obj2)
    if any(size(obj2)~=[obj1.n 1])
        eid = 'conZono:SizesNotEqual';
        msg = ['If input argument obj2 is a vector, can only add a zonotope in ',num2str(obj1.n),...
            ' dimensions to a vector with dimensions [',num2str(obj1.n),' 1]!'];
        throwAsCaller(MException(eid,msg))
    end
else
    throwAsCaller(MException('conZono:WrongInputType','Input argument obj2 must be either a Zono object or a column vector of size n.'));
end
end