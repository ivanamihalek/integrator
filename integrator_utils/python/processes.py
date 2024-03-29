
import multiprocessing
import os

########################################
def get_process_id():

	return os.getpid()


# don't know how to do this if there are no other_args,
# except by passing it an empty list in th place
########################################
def round_robin_pll(number_of_chunks, embarassingly_pllbl_fn, list, other_args):
	list_per_process = []
	for process in range(number_of_chunks):
		list_per_process.append([])

	for element in list:
		idx = list.index(element)
		list_per_process[idx%number_of_chunks].append(element)

	# run
	processes = []
	for ps in range (number_of_chunks):

		process = multiprocessing.Process(target=embarassingly_pllbl_fn, args=(list_per_process[ps], other_args))
		try:
			process.start()
			processes.append(process)
		except:
			print("Error: unable to start process")
			return False

	return processes

########################################
def parallelize (number_of_chunks, embarassingly_pllbl_fn, list, other_args):


	if (number_of_chunks < 1):
		print("number of processs is expected to be >= 1")
		return False

	if (number_of_chunks == 1):
		if other_args==None:
			ret = embarassingly_pllbl_fn (list)
		else:
			ret = embarassingly_pllbl_fn (list, other_args)
		return ret


	#########################################
	# nontrivial
	load = []
	for thr in range(number_of_chunks):
		load.append(0)

	for job in range (len(list)):
		load[(job%number_of_chunks)] += 1

	# run
	total = 0
	for ps in range (number_of_chunks):
		ps_from = total
		ps_to   = total + load[ps]
		total   += load[ps]

		if (ps_from >= len(list)):
			break
		if (ps == number_of_chunks-1):
			ps_to = len(list)


		process = multiprocessing.Process (target=embarassingly_pllbl_fn, args=(list[ps_from:ps_to], other_args))
		try:
			process.start()
		except:
			print("Error: unable to start process")
			return False
    
    
        
        
