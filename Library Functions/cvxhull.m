function [z] = cvxhull(x,y)
   % cvxhull computes the convex hull of two constrained zonotopes 
   % Inputs: x,y constrained zonotopes in CG-Rep as a struct variable with
   % c,G,A,b (c - center, G - generator)
   % satisfying x or y = {c + G\xi, ||\xi||_{\infty} \leq 1, A\xi = b (constraints)}.

   % Returns a constrained zonotope(z) in CG-Rep as a struct variable with
   % c,G,A,b

    n_dim = size(x.G,1); % Dimension
    n_c_x = size(x.A,1); % No. of constraints in x
    n_g_x = size(x.G,2); % No. of generators in x
    n_c_y = size(y.A,1); % No. of constraints in y
    n_g_y = size(y.G,2); % No. of generators in y
    
    z.c = (x.c + y.c)/2; % Computes the new center
  
    z.G = [x.G y.G (x.c-y.c)/2 zeros(n_dim,2*(n_g_x+n_g_y))]; % Computes the generators of z.
    z.H = [eye(n_g_x) zeros(n_g_x,n_g_y) -0.5*ones(n_g_x,1);... % Equality constraints
           -eye(n_g_x) zeros(n_g_x,n_g_y) -0.5*ones(n_g_x,1);...
           zeros(n_g_y,n_g_x) eye(n_g_y)  0.5*ones(n_g_y,1);...
           zeros(n_g_y,n_g_x) -eye(n_g_y)  0.5*ones(n_g_y,1)];
    k_1 = 1; k_2 = 1; % Handles the case when x or y is a singleton.
    if (n_g_x == 0)
        k_1 = 0;
    elseif (n_g_y == 0)
        k_2 = 0; 
    end
        
    z.I = 1*eye(2*(k_1*n_g_x + k_2*n_g_y)); % Adds constraints
    z.f = -0.5*ones(2*(k_1*n_g_x + k_2*n_g_y),1); % Satisfies [z.H z.I] \xi = z.f        
    z.A = [x.A zeros(n_c_x,n_g_y) -x.b/2 zeros(n_c_x,size(z.G,2)-n_g_x-n_g_y-1);... % Includes original constraints of x,y along with additional constraints added in Lines 20-23, 31-32.
           zeros(n_c_y,n_g_x) y.A +y.b/2 zeros(n_c_y,size(z.G,2)-n_g_y -n_g_x-1);...
           z.H z.I];
    z.b = [x.b/2;y.b/2;z.f];            
        
end