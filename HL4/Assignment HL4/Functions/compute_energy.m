function [energy] = compute_energy(input)   %function to compute the energy of a discrte-time signal
    len = length(input);
    energy = 0;

    for i = 1:len
        energy = energy + abs(input(i))^2;
    end

end

