N = 4;
BS = 2;

A = rand(N,N);
Acopy = A;
NB = N/BS;

mydet = eye(BS);

z = zeros(BS,BS);
one = eye(BS);

for k = 1:NB
   x = (k-1)*BS + 1;
   y = k*BS;
   pivot = A(x:y,x:y);
   mydet = pivot * mydet;
   A(:,x:y) = -[A(1:x-1,x:y);z;A(y+1:N,x:y)] / pivot
   A      = A + A(:,x:y)*A(x:y,:)
   A(x:y,:) = ([A(x:y,1:x-1),one,A(x:y,y+1:N)]' /pivot')'
end

A*Acopy
det(mydet)
det(Acopy)
