function i_bw = pectoral(w, h, a, b)
    i_bw = poly2mask([0 0 (-1 * b / a)], [0 b 0], w, h);
end
