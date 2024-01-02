function printToLatexFile(objectsArray, filename, xsol)
    fid = fopen(filename, 'w');
    
    if fid == -1
        error('No file.');
    end
    fprintf(fid, 'Name & Method & ParName & $L$ & $p$ & $i^*_{\\mathrm{L curve}}$ & $i^*_{\\mathrm{GCV}}$ & $\\left\\|\\mathbf{x}^* - \\mathbf{x}_{\\mathrm{L curve}}\\right\\|$ & $\\left\\|\\mathbf{x}^* - \\mathbf{x}_{\\mathrm{GCV}}\\right\\|$ \\\\\n');
    fprintf(fid, '\\hline\n');
    
    for i = 1:numel(objectsArray)
        name = objectsArray(i).name;
        method = objectsArray(i).method;
        parname = objectsArray(i).parname;
        Lname = objectsArray(i).Lname;
        ind_best_L = objectsArray(i).ind_best_L;
        ind_best_gcv = objectsArray(i).ind_best_gcv;
        p = max(objectsArray(i).p);
        x_L = objectsArray(i).Xs(:,objectsArray(i).ind_best_L);
        x_GCV = objectsArray(i).Xs(:,objectsArray(i).ind_best_gcv);
        fprintf(fid, '\n%s & \\texttt{%s} & $%s$ & $%s$ & %i & %d & %d & %d& %d\\\\\\hline',...
            name, method, parname, Lname, p, ind_best_L, ind_best_gcv, norm(x_L-xsol), norm(x_GCV-xsol));
    end
    fclose(fid);
end