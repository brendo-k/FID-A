files = dir("../metabolites/");
for i = 1:length(files)
    if ~isequal(files(i).name, '.') && ~isequal(files(i).name, '..')
        val = load(string(files(i).folder) + "/" + string(files(i).name));
        fields = fieldnames(val);
        metabolites{i-2} = val.(fields{1});
    end
end

t1 = zeros(length(metabolites), 1);
t2 = zeros(length(metabolites), 1);
for i = 1:length(metabolites)
    f = @() sim_Hamiltonian(metabolites{i}, 3);
    f2 = @() sim_Hamiltonian2(metabolites{i}, 3);
    t1(i) = timeit(f);
    t2(i) = timeit(f2);
end

fprintf("average time for sim_Hamiltonian: %d\n", mean(t1))
fprintf("average time for sim_Hamiltonian2: %d\n", mean(t2))
