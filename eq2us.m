function F = eq2us(x,a,f,N,A,epsilon)
xeq1=-N(1)/2:N(1)/2-1;
xeq2=-N(2)/2:N(2)/2-1;
mu1=-log(epsilon)/(2*N(1)^2)+3/(4*N(1)*N(1))*log(A);
mu2=-log(epsilon)/(2*N(2)^2)+3/(4*N(2)*N(2))*log(A);

Te1=1/pi*sqrt(-mu1*log(epsilon)+(mu1*N(1))^2/4+log(A)^2/(4*N(1)^2));
Te2=1/pi*sqrt(-mu2*log(epsilon)+(mu2*N(2))^2/4+log(A)^2/(4*N(2)^2));
M1=ceil(2*N(1)*Te1);
M2=ceil(2*N(2)*Te2);
phi0=zeros(2*N(1),2*N(2));
for i2=1:N(2)
    for i1=1:N(1)
        phi0(N(1)/2+i1,N(2)/2+i2)=exp(-mu1*xeq1(i1)^2-mu2*xeq2(i2)^2);
    end
end

fe=zeros(2*N(1),2*N(2));fe(N(1)/2+(1:N(1)),N(2)/2+(1:N(2)))=f./(2*N(1)*2*N(2))./phi0(N(1)/2+(1:N(1)),N(2)/2+(1:N(2)));
Fe=fftshift(fft2(fftshift(fe)));

%wrap
[idx,idy]=ndgrid(-M1:2*N(1)+M1-1,-M2:2*N(2)+M2-1);
idx0=mod(idx+2*N(1),2*N(1))+1;idy0=mod(idy+2*N(2),2*N(2))+1;
Fe(idx+M1+1,idy+M2+1)=Fe(idx0,idy0);

for k=1:size(x,1)
  F(k,1)=0;
  ell1=(floor(2*N(1)*x(k,1))-M1:floor(2*N(1)*x(k,1))+M1)';
  ell2=(floor(2*N(2)*x(k,2))-M2:floor(2*N(2)*x(k,2))+M2)';
  for i2=1:2*M2+1
    for i1=1:2*M1+1
      F(k,1)=F(k,1)+Fe(N(1)+M1+1+ell1(i1),N(2)+M2+1+ell2(i2)).*pi/(sqrt(mu1*mu2))*(exp(-pi^2/mu1*(ell1(i1)/(2*N(1))-x(k,1)-i*a(k,1)/(2*pi)).^2+ ...
                                                               -pi^2/mu2*(ell2(i2)/(2*N(2))-x(k,2)-i*a(k,2)/(2*pi)).^2));
    end
  end
end
