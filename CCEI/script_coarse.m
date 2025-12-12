addpath('../data');

filename = 'Halevy et al (2016) - Data';

subset = 3;

% import data
raw_data = readtable(strcat(filename, '.csv'));

raw_data.p_x = 1./table2array(raw_data(:, 'X_intercept'));
raw_data.p_y = 1./table2array(raw_data(:, 'Y_intercept'));

raw_data.prop_x = raw_data.X/(raw_data.X+raw_data.Y);

odds = raw_data.prop_x./(1-raw_data.prop_x);

id_vector = unique(raw_data.Subject)

n_subject = size(id_vector, 1);

results = zeros(n_subject, 3);

results(:, 1) = id_vector;

for i=1:n_subject
    i;
    id = id_vector(i)
    prices = raw_data(raw_data.Subject==id, {'p_x', 'p_y'});
    choices = raw_data(raw_data.Subject==id, {'X', 'Y'});

    prices = prices(1:subset, :);
    choices = choices(1:subset, :);

    results(i, 2) = CCEI(table2array(choices), table2array(prices), 0);
    results(i, 3) = CCEI(table2array(choices), table2array(prices), 1);
end

profiles.e_star = results(:, 2);
profiles.e_double_star = results(:, 3);

min(results(:, 2))

% writematrix(results, 'coarse_consistency_scores.csv')