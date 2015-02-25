function [vector]=mat_2_vect(Matrix)
[szm1, szm2]=size(Matrix);
inc_loc=1;
for j=1:szm2
    for i=1:szm1
        vector(inc_loc)=Matrix(i,j);
        inc_loc=inc_loc+1;
    end
end