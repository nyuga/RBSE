function [xrec,D] = initBinarySpectral(b,A)
if ~isreal(A)
    A=real(A);
end
y=b.^2;
Nin=size(A,2);
Nout=size(A,1);
p=mean((A),'all');

Ap=(1/p)*A-1;
Ea2=-1+1/p;
Ea3=1/(p^2)-3/p+2;
Ea4=1/(p^3)-4/(p^2)+6/p-3;
% Ea3=mean(Ap.^3,'all');
% Ea4=mean(Ap.^4,'all');
B=Ap'*((y/Nout).*(Ap))/(p^2);
% B=Ap'*diag(y/Nout)*(Ap)/(p^2);

B=B.*(1-eye(Nin));
C=mean(y.*Ap,1)/(p^2);
C=C(:);
D=B/Ea2^2 + 2/(Ea3+2*Ea2)*(diag(C-sum(B,2)/Ea2));

[eigDvec,eigDval]=eig(D);
[~,ind] = sort(diag(eigDval),'descend');
eigDvec=eigDvec(:,ind);
eigDvals=diag(eigDval(ind,ind));
% [eigDvec,eigDval]=eigs(D,2);
xrec=eigDvec(:,1)+1j*eigDvec(:,2);

xTx=(eigDvals(1)-eigDvals(2))/2;
xHx=(eigDvals(1)+eigDvals(2))/2;
% if xTx/xHx >0.1
% xrec=sqrt(xTx*xHx+xTx^2)*eigDvec(:,1)+sqrt(xTx*xHx-xTx^2)*1j*eigDvec(:,2);
% end
% xrec=xrec/norm(xrec)*sqrt(xHx);

Axrec=abs(A*xrec);
s = Axrec'*b/(Axrec'*Axrec);
xrec=xrec/norm(xrec)*s;
end
