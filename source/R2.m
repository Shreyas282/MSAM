function out = R2(x,y)

t = cov(x,y);
out = t(1,2)^2/(var(x)*var(y));
end