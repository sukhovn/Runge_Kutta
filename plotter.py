import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_data_file(file_name, legend=True):
	data = pd.read_csv(file_name)
	var_name = data.columns[0]
	func_names = data.columns[1:]

	var = data[var_name].to_list()
	funcs = []
	for name in func_names:
		funcs.append(data[name].to_list())

	plt.xlabel(var_name)
	for func in funcs:
		plt.plot(var, func)

	if(legend):
		plt.legend(func_names, loc="upper left")
	plt.show()

plot_data_file("Data/solution.dat")
# plot_data_file("Data/accuracy_st.dat")
# plot_data_file("Data/accuracy_int.dat")