if ~exist("sim_Hamiltonian" , 'file') || ~exist("sim_Hamiltonian2", 'file')
    error("please add FID-A and all subdierctories to your path")
end

result = runtests('compare_ham');
table(result)