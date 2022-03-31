function y = sem(x)
    y = std(x,'omitnan')./sqrt(nnz(~(isnan(x))));
end