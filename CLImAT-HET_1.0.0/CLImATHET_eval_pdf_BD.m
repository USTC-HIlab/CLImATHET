function results = CLImATHET_eval_pdf_BD(Bcount,Tcount,p)

if size(Bcount,1)>size(Bcount,2) %Nx1->1xN
    Bcount = Bcount';
end

if size(Tcount,1)>size(Tcount,2) %Nx1->1xN
    Tcount = Tcount';
end

if size(p,1)>size(p,2) %Nx1->1xN
    p = p';
end

results = binopdf(Bcount,Tcount,p);