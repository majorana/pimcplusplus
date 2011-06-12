x = rand(100,1);
N = length(x);

numk = 20;

Xc = zeros(numk,1);
Xs = zeros(numk,1);

t = [0:N-1]/N;

for ki=[0:numk-1]
  Xc(ki+1) = 0;
  Xs(ki+1) = 0;
  for i=[0:N-1]
%    Xc(ki+1) = Xc(ki+1) + 2.0*(x(i)*cos(2.0*pi*ki*t(i))/N);
%    Xs(ki+1) = Xc(ki+1) + 2.0*(x(i)*sin(2.0*pi*ki*t(i))/N);
     Xc(ki+1) = Xc(ki+1) + x(i+1)*exp(-sqrt(-1)*2.0*pi*ki*i/N)/N;
     
  end;
end;
%  Xc(1) = 0.5*Xc(1);
%  Xs(1) = 0;    
  
t2 = [0:999]/1000;
x2 = zeros(1000,1);

for ti =[1:1000]
  for ki=[0:numk-1]
%    x2(ti) = x2(ti) + Xc(ki+1)*cos(-2.0*pi*ki*t2(ti)); +...
%             Xs(ki+1)*sin(2.0*pi*ki*t2(ti));
    x2(ti) = x2(ti) + Xc(ki+1)*exp(sqrt(-1)*2.0*pi*ki*t2(ti));
  end;
end;

t3 = [0:1999]/1000;
x3 = [x2 ;x2];

plot (t, x, t3, real(x3));
