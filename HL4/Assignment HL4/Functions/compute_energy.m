function [energy] = compute_energy(input)
len = length(input);
energy = 0;

for i = 1:len
    energy = energy + abs(input(i))^2;
end

end

