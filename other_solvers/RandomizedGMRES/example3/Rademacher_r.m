function theta = Rademacher_r(k,n)
theta = 2/sqrt(k)*(rand(k,n)>1/2)-1/sqrt(k);
end
