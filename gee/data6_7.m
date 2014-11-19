data = csvread("epilepsy");
base_age_raw = dlmread ("epilepsy_base_age.txt", sep = " ");
base = real(base_age_raw);
age = imag(base_age_raw);
% repeat four times