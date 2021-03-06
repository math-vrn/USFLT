function F = eq2us_formula(x,a,f,N)
xeq1=-N(1)/2:N(1)/2-1;
xeq2=-N(2)/2:N(2)/2-1;
F=zeros(size(x,1),1);
for k=1:size(x,1)
    for i2=1:N(2)
        for i1=1:N(1)
          F(k) = F(k)+f(i1,i2)*exp(-2*pi*i*(xeq1(i1)*x(k,1)+xeq2(i2)*x(k,2))+...
                                        a(k,1)*xeq1(i1)+a(k,2)*xeq2(i2));
        end
    end
end