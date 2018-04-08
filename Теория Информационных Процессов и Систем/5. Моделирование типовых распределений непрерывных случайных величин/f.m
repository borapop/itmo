function p = f(x)
    p = 0;
    if sym(x <= 2)
        p = 0.3 .* x;
    end
    if sym(x > 2)
        p = 2;
    end
end