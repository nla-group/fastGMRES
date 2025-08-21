function m = getmem(w)
% return total memory in bytes used by all variables


m = 0;
for i = 1:length(w)
    m = m + w(i).bytes;
end

