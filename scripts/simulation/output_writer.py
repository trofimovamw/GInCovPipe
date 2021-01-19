
import pandas as pd

from pathlib import Path

import datetime
strptime = datetime.datetime.strptime

class writer:
	'''
	Class containing input and output information, reading and writing functions.
	'''

	def __init__(self, outputpath, file_prefix, init_date="2020-01-01"):
		self.fasta_path = outputpath + "/fasta"
		self.table_path = outputpath + "/table"
		self.init_date = init_date
		self.file_prefix = file_prefix

		try:
			Path(self.fasta_path).mkdir(parents=True, exist_ok=True)
		except FileExistsError:
			print("Creation of the directory %s failed" % self.fasta_path)
			exit()
		try:
			Path(self.table_path).mkdir(parents=True, exist_ok=True)
		except FileExistsError:
			print("Creation of the directory %s failed" % self.fasta_path)
			exit()


	def write_fasta(self, header_prefix, species_dict, sub_abs = None, file_suffix = None):
		'''
		Write the simulated sequences into one fasta file.
		Returns a table with dates and counts.
		'''
		outputfile = self.fasta_path + "/" + self.file_prefix + file_suffix + ".fasta"
		rows_list = []

		start_time = strptime(self.init_date, "%Y-%m-%d").date()
		t = 0
		try:
			file = open(outputfile, "w+")

			print("--- Write sequences into file " + outputfile + "---")

			for spec_dict in species_dict:
				date = start_time + datetime.timedelta(days=t)

				header = header_prefix + str(date.strftime("%Y-%m-%d"))

				sampled_N = 0
				for seq, num in spec_dict.items():
					if sub_abs is not None:
						# take abolsute subsample or N(t) if less
						num = min(sub_abs, num)
					sampled_N += num
					# print( (header+"\n"+seq+"\n") * num )
					file.write((header + "\n" + seq + "\n") * num)

				# faster way to add rows to a df
				dict1 = {"t": t, "date": date, "sampled_N": sampled_N}
				rows_list.append(dict1)

				t += 1

		except OSError:
			print("Writing of fasta file %s failed" % outputfile)
			exit()
		return pd.DataFrame(rows_list)

	def write_table(self, table, file_suffix = None):
	    '''
		Write the table with true and sampled sequence counts for each time step
		'''
	    outputfile = self.table_path + "/" + self.file_prefix + file_suffix + ".tsv"
	    print("--- Write table into file " + outputfile + "---")
	    table.to_csv(outputfile, sep="\t", header=True, index=False)
