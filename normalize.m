function J = normalize(I)
    J = I;
    minJ = min(J(:));
    maxJ = max(J(:));
    if maxJ > minJ
        J = (J-minJ)/(maxJ-minJ);
    end
end