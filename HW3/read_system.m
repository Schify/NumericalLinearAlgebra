function [arg1, arg2, arg3, arg4, arg5] = read_system(which_system)
% reads in one of the systems from the folders
% Which_system : 1 = Spiral Inductor
%                       [A, B, E]
%                2 = Butterfly Gyroscope
%                       [B, C, K, M, C_names]
%                3 = International space station component
%                       [A, B, C, D]
switch which_system
    case 1 % Spiral Inductor
        arg1 = mmread("Systems/PEECmodel/spiral_inductor_peec.A");
        arg2 = mmread("Systems/PEECmodel/spiral_inductor_peec.B");
        arg3 = mmread("Systems/PEECmodel/spiral_inductor_peec.E");

    case 2 % Butterfly Gyroscope
        arg1 = mmread("Systems/ButterflyGyro-dim1e5-gyro/gyro.B");
        arg2 = mmread("Systems/ButterflyGyro-dim1e5-gyro/gyro.C");
        arg3 = mmread("Systems/ButterflyGyro-dim1e5-gyro/gyro.K");
        arg4 = mmread("Systems/ButterflyGyro-dim1e5-gyro/gyro.M");
        arg5 = readtable("Systems/ButterflyGyro-dim1e5-gyro/gyro.C.names", ...
            "FileType", "text","Delimiter","\n","TextType","string","ReadVariableNames",false).Var1;

    case 3 % International space station component
        all_vars = load("Systems/Iss12a/iss12a.mat");
        arg1 = all_vars.A;
        arg2 = all_vars.B;
        arg3 = all_vars.C;
        arg4 = all_vars.D;
end



end

