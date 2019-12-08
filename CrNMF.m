function [B,H,F]=NNMF(W,Iter,k,beta,gamma)
% --Input  
%   --W is a [n,n,T] matrix, n is node number, T is the total time
%   --k represents dimension we select default 100
%   --beta is a parameter  default 1
%   --gama is a parameter default 1
%   --Iter is a parameter represents the iteration times default 300
%   author: Xiaoke Ma (xkma@xidian.edu.cn)



% --Output
%   --B is a [n,k,T] matrix
%   --H is a [k,k,T] matrix 
%   --F is a [k,n,T] matrix 

if nargin>5
	error('parameter is too much,the max number of parameter is 5');
end
switch nargin
	case 1
		Iter=300;
		k=300;
		beta=1;
		gamma=1;
	case 2
		k=300;
		beta=0.2;
		gamma=1;
	case 3 
		beta=1;
		gamma=1;
	case 4
		gamma=0.2;
end

[n,~,T]=size(W);
We=zeros(n,n,T);
%Get the degree matrix D
D=zeros(n,n,T);
for i=1:T
D(:,:,i)=diag(sum(W(:,:,i),2));
end 

%Get the degree matrix De
De=zeros(n,n,T);



%Initialization B H F with SVD method 
B=zeros(n,k,T);
H=zeros(k,k,T);
F=zeros(k,n,T);



for o=1:T
	[b,h,f]=svds(W(:,:,o),k);
	B(:,:,o)=abs(b);
	H(:,:,o)=abs(h);
	F(:,:,o)=abs(f');
end

A=B;
G=H;
E=F;


%Updata Rules

for i=1:T
	if i==1		
			B(:,:,i)=A(:,:,i).*((W(:,:,i)*E(:,:,i)'*G(:,:,i)') ./ (A(:,:,i)*G(:,:,i)*E(:,:,i)*E(:,:,i)'*G(:,:,i)'+eps));
			F(:,:,i)=E(:,:,i).*((G(:,:,i)'*A(:,:,i)'*W(:,:,i)) ./ (G(:,:,i)'*A(:,:,i)'*A(:,:,i)*G(:,:,i)*E(:,:,i)+eps));
			H(:,:,i)=G(:,:,i).*((A(:,:,i)'*W(:,:,i)*E(:,:,i)') ./ (A(:,:,i)'*A(:,:,i)*G(:,:,i)*E(:,:,i)*E(:,:,i)'+eps));
            M = abs(corr(B(:,:,i)'));
            We(:,:,i) = M;
            De(:,:,i) = diag(sum(We(:,:,i),2));
    else
        for o=1:Iter 
			B(:,:,i)=A(:,:,i).*((W(:,:,i)*E(:,:,i)'*G(:,:,i)'+beta*W(:,:,i-1)*A(:,:,i)) ./ (A(:,:,i)*G(:,:,i)*E(:,:,i)*E(:,:,i)'*G(:,:,i)'+beta*D(:,:,i-1)*A(:,:,i)+eps));
			F(:,:,i)=E(:,:,i).*((G(:,:,i)'*A(:,:,i)'*W(:,:,i)+gamma*E(:,:,i)*We(:,:,i-1)) ./ (G(:,:,i)'*A(:,:,i)'*A(:,:,i)*G(:,:,i)*E(:,:,i)+gamma*E(:,:,i)*De(:,:,i-1)+eps));
			H(:,:,i)=G(:,:,i).*((A(:,:,i)'*W(:,:,i)*E(:,:,i)') ./ (A(:,:,i)'*A(:,:,i)*G(:,:,i)*E(:,:,i)*E(:,:,i)'+eps));
			if norm(W(:,:,i)-B(:,:,i)*H(:,:,i)*F(:,:,i),'fro')<0.01
				break
			else
				A(:,:,i)=B(:,:,i);
				G(:,:,i)=H(:,:,i);
				E(:,:,i)=F(:,:,i);
			end
        end
        %Using the current time base matrix to constract the We 
        M = abs(corr(B(:,:,i)'));
        We(:,:,i) = M;
        De(:,:,i) = diag(sum(We(:,:,i),2));
	end
end
end
