import numpy as np
import os
import sys
import math

fn = sys.argv[1]
fwn = sys.argv[2]

gnList= []

f = open('gene_names.txt')

for line in f:

	gnList.append(line.strip())

f.close()

f = open(fn)

dataDic = {}


for line in f:

	tL = line.split()

	gn = tL[0]
	chrom = tL[1]
	pos = tL[2]
	ref = tL[3]
	var = tL[4]
	tcgaid = tL[5]


	ccfList = tL[6:]

	if not tcgaid in dataDic.keys():
		dataDic[tcgaid] = {}

#	dataDic[tcgaid][gn] = [chrom, pos, snp_class, var_type, ref, var, alt_count, ref_count, ccf_max, ccfList]
	dataDic[tcgaid][gn] = [chrom, pos, ref, var, ccfList]

f.close()

fw = open(fwn, 'w')

for x1 in range(len(gnList)):

	for x2 in range(len(gnList)):

		if x1 == x2:
			continue

		dist_ccf_Dic = {}

		g1 = gnList[x1]
		g2 = gnList[x2]


		for tcgaid in dataDic.keys():

			if g1 in dataDic[tcgaid].keys() and g2 in dataDic[tcgaid].keys():


				g1_ccf = dataDic[tcgaid][g1][-1]
				g2_ccf = dataDic[tcgaid][g2][-1]

				dist_ccf = []
				for x in range(199):
					dist_ccf.append(0.0)

				for x in range(-99, 100):

					idx = x
					sum_prob = 0.0

					for k in range(0, 100):
						g2_idx = int(-x+k)

						if g2_idx >= 0 and len(g2_ccf) > g2_idx:

							prob_g1 = g1_ccf[k]
							prob_g2 = g2_ccf[g2_idx]

							sum_prob += float(prob_g1) * float(prob_g2)

				
					dist_ccf[idx+99] = sum_prob


				max_dist_ccf_idx = 0

				for x in range(1, 199):
					if dist_ccf[x] > dist_ccf[max_dist_ccf_idx]:
						max_dist_ccf_idx = x

				max_dist_ccf_idx_null = 99
				for y in range(99, 199):
					if dist_ccf[y] > dist_ccf[max_dist_ccf_idx_null]:
						max_dist_ccf_idx_null = y

				if max_dist_ccf_idx < 99:
					dist_ccf_Dic[tcgaid] = [max_dist_ccf_idx, dist_ccf, max_dist_ccf_idx_null]





		loglik = 0.0

		for tcgaid in dist_ccf_Dic.keys():

			loglik += (math.log(dist_ccf_Dic[tcgaid][1][dist_ccf_Dic[tcgaid][2]]) - math.log(dist_ccf_Dic[tcgaid][1][dist_ccf_Dic[tcgaid][0]]))

		fw.write(g1 + "\t" + g2 + "\t" + str(loglik) + "\n")


fw.close()
