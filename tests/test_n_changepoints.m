clear
addpath(genpath('utils'));


tau_collect_test = [0 10 50 0 10 0 0 20 0  ];
[n_changepoints] = n_changepoints_from_list(tau_collect_test);

assert(all(n_changepoints==[2 1 0 1 0]))