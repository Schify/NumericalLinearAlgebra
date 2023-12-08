function printToLatexFile(objectsArray, filename, xsol)
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Check if the file was opened successfully
    if fid == -1
        error('Error opening the file.');
    end
    
    % Write LaTeX tabular header
%     fprintf(fid, '\\begin{tabular}{|c|c|c|c|c|}\n');
%     fprintf(fid, '\\hline\n');
    fprintf(fid, 'Name & Method & Parameter Name & $L$ & $p$ & $i^*_{\\mathrm{L curve}}$ & $i^*_{\\mathrm{GCV}}$ & $\\left\\|\\mathbf{x}^* - \\mathbf{x}_{\\mathrm{L curve}}\\right\\|$ & $\\left\\|\\mathbf{x}^* - \\mathbf{x}_{\\mathrm{GCV}}\\right\\|$ \\\\\n');
    fprintf(fid, '\\hline\n');
    
    % Iterate through each object in the array
    for i = 1:numel(objectsArray)
        % Extract properties from the object
        name = objectsArray(i).name;
        method = objectsArray(i).method;
        parname = objectsArray(i).parname;
        Lname = objectsArray(i).Lname;
        ind_best_L = objectsArray(i).ind_best_L;
        ind_best_gcv = objectsArray(i).ind_best_gcv;
        p = max(objectsArray(i).p);
        x_L = objectsArray(i).Xs(:,objectsArray(i).ind_best_L);
        x_GCV = objectsArray(i).Xs(:,objectsArray(i).ind_best_gcv);
        % Write the values to the file
        fprintf(fid, '\n%s & \\texttt{%s} & $%s$ & $%s$ & %i & %d & %d & %d& %d\\\\\\hline',...
            name, method, parname, Lname, p, ind_best_L, ind_best_gcv, norm(x_L-xsol), norm(x_GCV-xsol));
    end
    
    % Write LaTeX tabular footer
   % fprintf(fid, '\\hline\n');
    
    % Close the file
    fclose(fid);
end