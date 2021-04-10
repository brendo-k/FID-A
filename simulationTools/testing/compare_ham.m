clear
%preconditions
files = dir("../metabolites/");
for i = 1:length(files)
    if ~isequal(files(i).name, '.') && ~isequal(files(i).name, '..')
        val = load(string(files(i).folder) + "/" + string(files(i).name));
        fields = fieldnames(val);
        metabolites{i-2} = val.(fields{1});
    end
end

%% Test 1: Ix
for i = 1:length(metabolites)
    [H_one, d_one] = sim_Hamiltonian(metabolites{i}, 3);
    [H_two, d_two] = sim_Hamiltonian(metabolites{i}, 3);
    if length(H_one) > 1     
        for i = 1:length(H_one)
            assert(isequal(H_one(i).Ix , H_two(i).Ix))
        end
    else
        assert(isequal(H_one.Ix , H_two.Ix))
    end
end
%% Test 2: Ix
for i = 1:length(metabolites)
    [H_one, d_one] = sim_Hamiltonian(metabolites{i}, 3);
    [H_two, d_two] = sim_Hamiltonian(metabolites{i}, 3);
    if length(H_one) > 1     
        for i = 1:length(H_one)
            assert(isequal(H_one(i).Iy , H_two(i).Iy))
        end
    else
        assert(isequal(H_one.Iy , H_two.Iy))
    end
end
%% Test 3: Iz
for i = 1:length(metabolites)
    [H_one, d_one] = sim_Hamiltonian(metabolites{i}, 3);
    [H_two, d_two] = sim_Hamiltonian(metabolites{i}, 3);
    if length(H_one) > 1     
        for i = 1:length(H_one)
            assert(isequal(H_one(i).Iz , H_two(i).Iz))
        end
    else
        assert(isequal(H_one.Iz , H_two.Iz))
    end
end
%% Test 4: Fx
for i = 1:length(metabolites)
    [H_one, d_one] = sim_Hamiltonian(metabolites{i}, 3);
    [H_two, d_two] = sim_Hamiltonian(metabolites{i}, 3);
    if length(H_one) > 1     
        for i = 1:length(H_one)
            assert(isequal(H_one(i).Fx , H_two(i).Fx))
        end
    else
        assert(isequal(H_one.Fx , H_two.Fx))
    end
end
%% Test 5: Fy
for i = 1:length(metabolites)
    [H_one, d_one] = sim_Hamiltonian(metabolites{i}, 3);
    [H_two, d_two] = sim_Hamiltonian(metabolites{i}, 3);
    if length(H_one) > 1     
        for i = 1:length(H_one)
            assert(isequal(H_one(i).Fy , H_two(i).Fy))
        end
    else
        assert(isequal(H_one.Iz , H_two.Iz))
    end
end
%% Test 6: Fz
for i = 1:length(metabolites)
    [H_one, d_one] = sim_Hamiltonian(metabolites{i}, 3);
    [H_two, d_two] = sim_Hamiltonian(metabolites{i}, 3);
    if length(H_one) > 1     
        for i = 1:length(H_one)
            assert(isequal(H_one(i).Fz , H_two(i).Fz))
        end
    else
        assert(isequal(H_one.Iz , H_two.Iz))
    end
end
