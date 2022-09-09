function A = build_A(T,K,V,t,h)

L = length(h);
dh = h(2)-h(1);

lv = localvol(T,K,V,t,exp(h));
lvsq = lv.*lv;

a=zeros(L,1);
b=zeros(L,1);
c=zeros(L,1);

a(1:end-2,1) = 0.5*lvsq(2:end-1)*( -1/(2*dh) - 1/(dh^2) );
b(2:end-1,1) = 0.5*lvsq(2:end-1)*( 2/(dh^2) ); 
c(3:end,1)   = 0.5*lvsq(2:end-1)*( 1/(2*dh) - 1/(dh^2) ); 

% init conditions: a(end)=0; b(1)=0; b(end)=0; c(1)=0;

A = spdiags([a b c],[-1 0 1],L,L); 
%full(A);

end
