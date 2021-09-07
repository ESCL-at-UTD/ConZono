function out = volume(obj)


if obj.nC == 0
    % Adpted from CORA 2020 and from
    % [1] Filliman, Paul. "Extremum problems for zonotopes." Geometriae Dedicata 27.3 (1988): 251-262.
    if obj.nG<obj.n
        out=0;
        return;
    end
    
    if true
        % this does not work for very large number of
        % generators/dimensions
        % ==> runs out of memory. Although in these cases the number of
        % combinations is so large or computation of determinants so costly that it may take too long to compute
        % the volume. It may make more sense to compute a faster
        % approximation instead.
        indexes = nchoosek(1:obj.nG,obj.n);
        vol = 0;
        for i=1:size(indexes,1)
            A = obj.G(:,indexes(i,:));
            vol = vol + abs(det(A));
        end
        
    else % DO NOT USE THIS OPTION, RESULTS ARE SLIGHTlY OFF, NEEDS MORE DEBUGGING
        %This option does not run out of memory when computing the
        %permutations of generators since it computes one permutation
        %at a time.
        
        %auxiliar variables used to generate permuations of generators from
        %the generators matrix
        [a,indexes,x,y,z,~,p] = initializeTwiddle(obj.nG,obj.n);
        
        vol = 0;
        while true
            A = obj.G(:,indexes);
            vol = vol + abs(det(A));
            
            %gets next permutation
            [done,p,x,y,z] = permutation_twiddle(x,y,z,p);
            if done
                break;
            else
                indexes(z+1) = a(x+1);
            end
            
        end
        
    end
    
    out = 2^obj.n*vol;
    return;
    
else
    Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1),'He',[obj.A obj.b]);
end
P = plus(obj.c,affineMap(Box,obj.G));
out = P.volume;
end




%% Adapted from http://www.netlib.no/netlib/toms/382 ( Coded by Matthew Belmonte <mkb4@Cornell.edu>, 23 March 1996)
% See https://stackoverflow.com/a/127856
function [done,p,x,y,z] = permutation_twiddle(x,y,z,p)

j=0;
j=j+1;
while p(j+1)<=0
    j = j+1;
end

if p(j)==0
    for i=j-1:-1:2
        p(i+1)=-1;
    end
    p(j+1)=0;
    p(2) = 1;
    x=0;
    z=0;
    y=j-1;
    
else
    
    if j>1
        p(j)=0;
    end
    
    %L2:
    j = j+1;
    while p(j+1)>0
        j = j+1;
    end
    i=j-1;
    k=j-1;
    
    %L3:
    i=i+1;
    while p(i+1)==0
        p(i+1)=-1;
        i = i+1;
    end
    
    if p(i+1)==-1
        p(i+1)=p(k+1);
        z=p(k+1)-1;
        x=i-1;
        y=k-1;
        p(k+1) = -1;
    else
        
        if i==p(1)
            done=true;
            return;
        else
            
            z = p(i+1)-1;
            p(j+1) = p(i+1);
            p(i+1) = 0;
            x=j-1;
            y=i-1;
        end
        
    end
end
done = false;
end


%%
function [a,c,x,y,z,done,p] = initializeTwiddle(n,m)

a = 1:n;
c = a(n-m+1:end);

done=false;
x=1;
y=1;
z=1;

%initialize p
p = zeros(n+2,1);
p(1) = n+1;
p(n+2) = -2;
% p(2:(n-m+1))=0;
p(n-m+2:n+1)=1:m;

end