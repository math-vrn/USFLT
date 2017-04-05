A=100;epsilon=1e-12;
N=[64 32];
M=prod(N);
x=(rand(M,2)-0.5);
a=(rand(size(x))-0.5)*log(A)./repmat(N,M,1)*2;

%us2eq
f=rand(M,1)+i*rand(M,1);
F_formula=us2eq_formula(x,a,f,N);
F=us2eq(x,a,f,N,A,epsilon);
norm(F-F_formula,'inf')
%eq2us
f=rand(N)+i*rand(N);
G_formula=eq2us_formula(x,a,f,N);
G=eq2us(x,a,f,N,A,epsilon);
norm(G-G_formula,'inf')