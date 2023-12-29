function lot(N, R)
%lot picks lottery numbers. N numbers from 1 to R are chosen.
if (N <= R & R > 0)
    tot = randperm(R);
%  returns  a  row  vector  containing  a  random  permutation  of  the  integers  from  1  to n without 
%repeating %elements.
    x = tot(1:N);
    disp(['Your lottery numbers are: ', num2str(x)])
end