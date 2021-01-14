#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from random import choice, choices
import numpy as np
import csv
import os
import datetime
import random
import math
import copy

from pathlib import Path

strptime = datetime.datetime.strptime


class sequenceEvol:
	'''
	Each instance of this class represens one evolutionary trajectory run. Starting with a sequence set of equal sequences,
	growing and mutating over time.
	'''

	def __init__(self, outFile, length, p_mut, N_init, t_start, t_final,  p_repl=None, p_repl2=None, p_death=None, t_switch=None, init_seq=None):
		self.L = length
		self.p_repl = p_repl
		self.p_repl2 = p_repl2
		self.p_death = p_death
		self.p_mut = p_mut
		self.N = N_init
		self.t_start = t_start
		self.t_final = t_final
		self.t_switch = t_switch
		self.outfile = outFile
		self.init_seq = self._initSeqArr(init_seq)
		#self.init_seqs = self._initSeqSet(init_seq)
		self.initB = self._initSetB()
		#self.sim = sim
		#self.total_sim = total_sim
		#for s in self.init_seqs:
		#		print(s)
		print(self.init_seq)

	def _initSeqArr(self, init_seq=None):
		'''
		Create initial sequence as array of each nucleotide. If no sequence is given, create a random one.
		'''
		if init_seq is None:
			alphabet = 'ACGT'
			init_seq = ''.join(random.SystemRandom().choice(alphabet)for _ in range(self.L))
		return list(init_seq)

	# def _initSeqSet(self, init_seq=None):
	# 	'''
	# 	Create initial sequence set of size N. If no sequence is given, create a random one.
	# 	'''
	# 	if init_seq is None:
	# 		alphabet = 'ACGT'
	# 		init_seq = ''.join(random.SystemRandom().choice(alphabet)for _ in range(self.L))
	# 	return [init_seq for _ in range(self.N)]

	def _initSetB(self):
		'''
		TODO: was tut es?
		'''
		z = np.zeros(self.L)
		zset = []
		for i in range(self.N):
			zset.append((i, z))
		return zset

	def _mutateBase(self, base, initSeqBase):
		'''
		TODO: Kommentare entfernen
		Mutate individual sequences according
		to Kimura model, fixed transition probability
		'''
		# Mutate if the base didnt mutate yet
		#if base == initSeqBase:
		stay = 0
		alpha = 0.2
		beta = 0.2
		gamma = 0.6

		transitions = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
		transversions_alpha = {'A': 'C', 'C': 'A', 'T': 'G', 'G': 'T'}
		transversions_beta = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}

		r = np.random.uniform(0, 1)

		if stay <= r < (stay + gamma):
			return transitions[base]
		elif (stay + gamma) <= r < (stay + gamma + alpha):
			return transversions_alpha[base]
		elif (stay + gamma + alpha) <= r < (stay + gamma + alpha + beta):
			return transversions_beta[base]
		return base
		#else:
			#return base

	def mutate(self):
		'''
		Mutate initial sequence set over the course of t generations
		'''
		init = self.init
		print(init)
		initSeq = init[0]


		# start time
		t = self.t_start

		# collect samped time steps
		time_steps = []
		time_steps.append(init)

		trajectory = []
		time_trajectory = []
		tipsList = []

		# the current replication rate of the time step
		curr_p_repl = self.p_repl

		while t < self.t_final:
			print("Current time t = ", t)

			# Switch the replication rate somewhere in the middle if that is the correct mode (p_repl2!=0)
			if (self.p_repl2 is not None) and (t >= self.t_switch):
				curr_p_repl = self.p_repl2

			r_repl = curr_p_repl * len(init)
			# print("Replication rate: ",r_repl)
			r_death = self.p_death * len(init)
			# print("Death rate: ",r_death)
			rate = (r_repl, r_death)

			u1 = random.uniform(0, 1)
			u2 = random.uniform(0, 1)
			while isinstance(u1, complex):
				u1 = random.uniform(0, 1)
			delta_t = (1 / sum(rate)) * math.log(1 / u1)
			while isinstance(u2, complex):
				u2 = random.uniform(0, 1)

			csum_rate = (r_repl / (r_repl + r_death), r_death / (r_repl + r_death))
			csum = np.cumsum(csum_rate)
			r = u2 * sum(csum_rate)

			repl_event = False
			death_event = False

			if r < csum[0]:
				repl_event = True
			elif csum[0] < r < csum[1]:
				death_event = True

			new_set = []

			# Add sequences to list if birth event
			if repl_event == True:
				# print("  Picked replication event")
				# Pick sequences to replicate with a certain probability - by index
				pick_ind = random.choices(np.arange(0, len(init)), weights=np.full(len(init), curr_p_repl / r_repl),
				                          k=math.ceil(len(init) * (curr_p_repl / r_repl)))
				# print("Sequences to replicate: ",pick_ind)
				# Initialize mutated sites
				num_mut_sites = len(pick_ind) * L + 1
				while num_mut_sites >= len(pick_ind) * L:
					num_mut_sites = np.random.poisson(self.p_mut * len(pick_ind) * L)
				# print("                Number of mutated sites: ",num_mut_sites)
				ind_mut_sites = np.arange(0, len(pick_ind) * L)
				mut_sites = random.sample(ind_mut_sites.tolist(), num_mut_sites)

				pick_children = []
				for ind in pick_ind:
					pick_children.append(init[ind])
					pick_children.append(init[ind])

				if num_mut_sites > 0:
					for i in mut_sites:
						seq_ind = math.floor(i / L)
						seq_pos = i % L
						seq = list(pick_children[seq_ind])
						# seq[seq_pos] = self._mutateBase(new_set[seq_ind][seq_pos],initSeq[seq_pos])
						seq[seq_pos] = self._mutateBasenorm(pick_children[seq_ind][seq_pos])
						# print("                       Mutated site: ",seq_pos)
						s = "".join(seq)
						pick_children[seq_ind] = s
				# If birth event - no tips available
				tipsList.append([])
				# Append children to init
				for seq in pick_children:
					init.append(seq)
				# print("   New init size: ",len(init))
			# if death event - remove sequences from population
			elif death_event == True:
				# print("Picked death event")
				pick = []
				pick_ind = random.choices(np.arange(0, len(init)), weights=np.full(len(init), self.p_death / r_death),
				                          k=math.ceil(len(init) * (self.p_death / r_death)))

				new_init = []
				for idx in range(len(init)):
					if idx in pick_ind:
						pick.append(init[idx])
					else:
						new_init.append(init[idx])
				init = new_init
				# print("   New init size: ",len(init))
				tipsList.append(pick)

			# Add to trajectory
			new_set = copy.copy(init)
			trajectory.append(new_set)
			t = t + delta_t
			time_trajectory.append(t)

			if not new_set:
				break

		return trajectory, time_trajectory, initSeq, tipsList
