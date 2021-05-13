import os

number_of_fits = 1000

fit_number_list_file = "fit_number_list.txt"

if os.path.exists(fit_number_list_file):
    os.remove(fit_number_list_file)

with open(fit_number_list_file, "w") as fit_list:
    for i in range(number_of_fits):
        if not os.path.exists('output/fit_result_{:06d}.root'.format(i)):
            fit_list.write('{:d} fit_result_{:06d}.root\n'.format(i, i))

