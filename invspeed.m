function lot(3, 6000)
%lot picks lottery numbers. N numbers from 1 to R are chosen.
if (3 <= 600 & 600 > 0)
    tot = randperm(600);
%  returns  a  row  vector  containing  a  random  permutation  of  the  integers  from  1  to n without 
repeating %elements.
    x = tot(1:3);
    disp(['Your lottery numbers are: ', num2str(x)])
end
