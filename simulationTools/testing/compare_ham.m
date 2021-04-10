clear
%preconditions
files = dir("../metabolites/");
for i = 1:length(files)
    if ~isequal(files(i).name, '.') && ~isequal(files(i).name, '..')
        metabolites{i-2} = load(string(files(i).folder) + "/" + string(files(i).name));
    end
end

%% Test 1: Ix
for i = 1:length(metabolites)
    [H_one, d_one] = sim_Hamiltonian(metabolites{i}, 3);
    [H_two, d_two] = sim_Hamiltonian(metabolites{i}, 3);
    assert(isequal(H_one.Ix , H_two.Ix))
end
