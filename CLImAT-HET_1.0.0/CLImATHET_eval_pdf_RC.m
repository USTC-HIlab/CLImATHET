function results = CLImATHET_eval_pdf_RC(data,lambda,p)

if size(data,1)>size(data,2) %Nx1->1xN
    data = data';
end

if size(lambda,1)>size(lambda,2) %Nx1->1xN
    lambda = lambda';
end

if size(p,1)>size(p,2) %Nx1->1xN
    p = p';
end
temp = gammaln(data+lambda.*(1-p)./p)+(lambda.*(1-p)./p).*log(1-p)+data.*log(p)...
    -gammaln(data+1)-gammaln(lambda.*(1-p)./p);
results = exp(temp);
results(results <=0) = eps;