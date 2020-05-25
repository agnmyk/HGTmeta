import numpy
import copy
import random
import itertools
import sys
import getopt
from tree import Node

def print_help():
	print "\nUsage:"
	print "hgtboot.py -S <file with a species tree> -G <file with a gene tree> [-D <int>, -L <int>, -T <int>, -h]"
	print "\nOptions:"
	print "\t-S <file with a species tree>"
	print "\t-G <file with a gene tree>"
	print "\t-D <cost of duplication event>  -- default value = 2"
	print "\t-L <cost of loss event> -- default value = 1"
	print "\t-T <cost of horizontal gene transfer event> -- default value = 3"
	print "\t-h -- print help"
	sys.exit()



Ploss = 1
Ptheta = 3
Pdelta = 2



try:
	opts, args = getopt.getopt(sys.argv[1:], "S:G:T:D:L:i:h")
except getopt.GetoptError as err:
	print "Wrong option"
	print_help()


for o, a in opts:
	if o == "-S": stree = a
	elif o == "-G": gtree = a
	elif o == "-D": Pdelta = int(a)
	elif o == "-T": Ptheta = int(a)
	elif o == "-L": Ploss = int(a)
	elif o == "-h": print_help()
	
if not any("-S" in a for a in opts): 
	print "No species tree given\n\n"
	print_help()

if not any("-G" in a for a in opts): 
	print "No gene tree given\n\n"
	print_help()

try:

	streeFin = open(stree)
	S = Node.readTree(streeFin.readline())
	gtreeFin = open(gtree)
	G = Node.readTree(gtreeFin.readline())

except:
	print "Could not load the gene tree or the species tree"
	sys.exit(2) 
	

size_S = len(S.leaves())+len(S.inner())
size_G= len(G.leaves())+len(G.inner())

def nodesByID(n, d):
	d[n.nr] = n
	if not n.leaf:
		for s in n.sons: nodesByID(s, d)
	return d

def find_all_mins(lst):
	mins = []
	minn = min(lst)
	prev_minn = minn
	c = 0
	while prev_minn == minn and c < len(lst):	
		ind = lst.index(minn)
		mins.append(ind)
		lst[ind] = float("inf")
		minn = min(lst)
		c += 1
	return mins
	
def count_results(X, Y):
	result = reconstruct[X][Y]
	if not result_counters[X][Y]:
		if result:
			for res in result:
				count_results(res[0][0], res[0][-1])
				count_results(res[1][0], res[1][-1])
			for res in result: 
				numOfCh = sum(result_counters[res[0][0]][res[0][-1]])* sum(result_counters[res[1][0]][res[1][-1]])
				
				result_counters[X][Y].append(numOfCh)
		else:	result_counters[X][Y].append(1)



def random_result(g, s):
	result = reconstruct[g][s]
	res = []
	if len(result) > 1:
		results = {k:[] for k in result_counters[g][s]}
		for i in xrange(len(result_counters[g][s])): 
			results[result_counters[g][s][i]].append(result[i])
		sorted_weights = sorted(results.keys())
		summ = sum(sorted_weights)
		i = random.random() * summ
		maxweight, wsum = 0, 0
		for w in sorted_weights:
			wsum += w
			if i < wsum:
				maxweight = w
				break
		i = random.randint(0, len(results[maxweight])-1)
		res = results[maxweight][i]
	elif len(result) == 1:
		res = result[0]
	nodeG = GnodesByID[g]
	nodeS = SnodesByID[s]
	
	if nodeG.leaf: 
		heatmap_results[nodeG.nr][nodeS.nr] += 1
	if len(result) >= 1:
		random_result(res[0][0], res[0][-1])
		random_result(res[1][0], res[1][-1])

SnodesByID = nodesByID(S, {})
GnodesByID = nodesByID(G, {})



G.lca(S)


leaves_G = G.leaves()
leaves_S = S.leaves()


niesk = float("inf")

c=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
c_sigma=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
c_delta=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
c_theta=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
inn=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
inn_rec = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]
unknownleaves_rec = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]
inAlt=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
inAlt_rec = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]
out=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
out_res = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]
events = numpy.array([['-']*(size_S+1)]*(size_G+1))
reconstruct = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]
result_counters = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]

heatmap_results = {k.nr:{kk.nr:0 for kk in S.leaves()+S.inner()} for k in G.leaves()}	

Result = [[]]
	

					

if True:	
	for g in leaves_G:
		if "?" in g.label:
			for i in xrange(len(c[g.nr])-1): 
				nodeS = SnodesByID[i]
				if nodeS.leaf: c[g.nr][i] = 0
		else: c[g.nr][g.map.nr] = 0

		if "?" not in g.label:
			inn[g.nr][g.map.nr] = Ploss*0 
			inn_rec[g.nr][g.map.nr] = [g.map.nr]
			inAlt[g.nr][g.map.nr] = 0
			inAlt_rec[g.nr][g.map.nr] = [g.map.nr]
			ps = g.map.fath
			dS = 1
			while ps != None:
				inn[g.nr][ps.nr] = Ploss*dS
				inn_rec[g.nr][ps.nr] = [g.map.nr]
				inAlt[g.nr][ps.nr] = 0
				inAlt_rec[g.nr][ps.nr] = [g.map.nr]
	    			ps=ps.fath
				dS += 1
		else:
			def findLeaf(snodes, depth):
				wasLeaf = False
				
				leavesp = []
				for snode in snodes:
					for snode_son in snode.sons:
						leavesp.append(snode_son)
						if snode_son.leaf: wasLeaf = True
				depth += 1	
				if wasLeaf:
					return [[l for l in leavesp if l.leaf], depth]
				else:
					return findLeaf(leavesp, depth)
				 

			for snode in S.inner():
				goodleaves, leavesdepth = findLeaf([snode], 0)
				inn[g.nr][snode.nr] = Ploss*leavesdepth
				unknownleaves_rec[g.nr][snode.nr] = [gl.nr for gl in goodleaves]
				inAlt[g.nr][snode.nr] = 0
				l =S.leaves()
				for sleaf in snode.leaves():
					l.remove(sleaf)
				inAlt_rec[g.nr][snode.nr] = [ll.nr for ll in l]

			for gmap in S.leaves():
				inn[g.nr][gmap.nr] = Ploss*0
				inAlt[g.nr][gmap.nr] = 0
				l =S.leaves()
				l.remove(gmap)
				inAlt_rec[g.nr][gmap.nr] = [ll.nr for ll in l]
				unknownleaves_rec[g.nr][gmap.nr] = [gmap.nr]

	g=G
	def cos(g, s):
		if not g.leaf:
			gp, gpp = g.sons
			cos(gp, s) #
			cos(gpp, s) #	
			def cos3(g, s):
				if not s.leaf:
					for sson in s.sons:
						other_ssons = [spson for spson in s.sons if spson.nr != sson.nr]
						costs = [out[g.nr][s.nr]] + [inAlt[g.nr][spson.nr] for spson in s.sons if spson.nr != sson.nr]
						out_g_sson = min(costs)
						inds = find_all_mins(costs)
						rec = []
						for ind in inds:
							if ind == 0:
								if s.fath:
									other_nodes = []
									fath = s.fath
									currentS = s
									
									while fath:
										other_nodes.extend([son for son in fath.sons if son.nr != currentS.nr])
										currentS = fath
										fath = fath.fath
									inds2 = find_all_mins([inAlt[g.nr][node.nr] for node in other_nodes])
									for ind2 in inds2:										
										rec.extend(inAlt_rec[g.nr][other_nodes[ind2].nr])
	
								else:
									continue
							else:
								rec.extend(inAlt_rec[g.nr][other_ssons[ind-1].nr])
								
						out[g.nr][sson.nr] = out_g_sson
						out_res[g.nr][sson.nr] = copy.deepcopy(rec)

					for sson in s.sons: cos3(g, sson)

			def cos2(g, s):		
				gfath = g.fath.nr if g.nr != 0 else 0
				sfath = g.fath.map.nr if g.nr != 0 else 0
				if s.leaf:
					c_sigma[g.nr][s.nr] = niesk
					c_delta[g.nr][s.nr] = Pdelta + c[gp.nr][s.nr] + c[gpp.nr][s.nr]
					reconstruct_delta = [((gp.nr,s.nr),(gpp.nr,s.nr), (gfath,sfath), "D")] ###
					if s.nr != 0: 
						c_theta[g.nr][s.nr] = Ptheta + min(inn[gp.nr][s.nr] + out[gpp.nr][s.nr], inn[gpp.nr][s.nr] + out[gp.nr][s.nr])
						inds = find_all_mins([Ptheta + inn[gp.nr][s.nr] + out[gpp.nr][s.nr],Ptheta + inn[gpp.nr][s.nr] + out[gp.nr][s.nr]])
						reconstruct_theta = []
						
						for ind in inds:
							if ind == 0:
								for out_result in out_res[gpp.nr][s.nr]:
									reconstruct_theta.append(((gp.nr,s.nr),(gpp.nr, "->", out_result), (gfath,sfath),"T")) 
							else:
								for out_result in out_res[gp.nr][s.nr]:
									reconstruct_theta.append(((gpp.nr,s.nr),(gp.nr, "->", out_result), (gfath,sfath),"T"))

					c[g.nr][s.nr] = min(c_delta[g.nr][s.nr], c_theta[g.nr][s.nr], c_sigma[g.nr][s.nr])
					inn[g.nr][s.nr] = c[g.nr][s.nr]
					inn_rec[g.nr][s.nr] = [s.nr]
					inAlt[g.nr][s.nr] = c[g.nr][s.nr]
					inAlt_rec[g.nr][s.nr] = [s.nr]

					inds = find_all_mins([c_delta[g.nr][s.nr], c_theta[g.nr][s.nr], c_sigma[g.nr][s.nr]])
					
					for ind in inds:
						if ind == 0: 
							reconstruct[g.nr][s.nr].extend(reconstruct_delta)
						elif ind == 1: 
							reconstruct[g.nr][s.nr].extend(reconstruct_theta)

				else:
					ssons_pairs = list(itertools.combinations(s.sons, 2))
					numOf_ssons = len(s.sons)
					for sson in s.sons: cos2(g, sson) 

					sigma_min = 100*niesk
					sigma_min_sons = []
					numOf_sigma_losses = numOf_ssons - 2 
					for sp, spp in ssons_pairs:
						minn = min(inn[gp.nr][sp.nr] + inn[gpp.nr][spp.nr], inn[gpp.nr][sp.nr] + inn[gp.nr][spp.nr])
						if minn < sigma_min: 
							sigma_min = minn
							sigma_min_sons = [(sp, spp)]
						elif minn == sigma_min: 
							sigma_min_sons.append((sp, spp))
					reconstruct_sigma = []
					c_sigma[g.nr][s.nr] = sigma_min + Ploss * numOf_sigma_losses
					for sp, spp in sigma_min_sons:

						inds = find_all_mins([inn[gp.nr][sp.nr] + inn[gpp.nr][spp.nr], inn[gpp.nr][sp.nr] + inn[gp.nr][spp.nr]])
						
						for ind in inds:
							if ind == 0: 
								gpmaplist = unknownleaves_rec[gp.nr][sp.nr]  if "?" in gp.label else inn_rec[gp.nr][sp.nr]
								gppmaplist = unknownleaves_rec[gpp.nr][spp.nr]  if "?" in gpp.label else inn_rec[gpp.nr][spp.nr]
								for mapspair in list(itertools.product(gpmaplist, gppmaplist)):
									reconstruct_sigma.append(((gp.nr,mapspair[0]),(gpp.nr,mapspair[1]), (gfath,sfath),"S"))
								
							else: 
								
								gpmaplist = unknownleaves_rec[gp.nr][spp.nr]  if "?" in gp.label else inn_rec[gp.nr][spp.nr]
								gppmaplist = unknownleaves_rec[gpp.nr][sp.nr]  if "?" in gpp.label else inn_rec[gpp.nr][sp.nr]
								for mapspair in list(itertools.product(gppmaplist, gpmaplist)):
									reconstruct_sigma.append(((gpp.nr,mapspair[0]),(gp.nr,mapspair[1]), (gfath,sfath),"S"))

					delta_min = 100*niesk
					delta_min_sons = []
					for sp, spp in ssons_pairs:
						for sp2, spp2 in ssons_pairs:
							minn = min(
										c[gp.nr][s.nr]+inn[gpp.nr][spp2.nr]+(numOf_ssons-1)*Ploss,
										c[gp.nr][s.nr]+inn[gpp.nr][sp2.nr]+(numOf_ssons-1)*Ploss,
										c[gpp.nr][s.nr]+inn[gp.nr][spp.nr]+(numOf_ssons-1)*Ploss,
										c[gpp.nr][s.nr]+inn[gp.nr][sp.nr]+(numOf_ssons-1)*Ploss,
										c[gp.nr][s.nr]+c[gpp.nr][s.nr],
										inn[gp.nr][sp.nr]+inn[gpp.nr][spp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss,
										inn[gp.nr][spp.nr]+inn[gpp.nr][sp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss,
										inn[gp.nr][sp.nr]+inn[gpp.nr][sp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss,
										inn[gp.nr][spp.nr]+inn[gpp.nr][spp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss)
							

							if minn < delta_min:
								delta_min_sons = [((sp, spp), (sp2, spp2))]
								delta_min = minn
								break
							
							elif minn == delta_min:
								delta_min_sons.append(((sp, spp), (sp2, spp2)))

					reconstruct_delta = []	

					c_delta[g.nr][s.nr] = Pdelta + delta_min
					added_pairs = {}
					for sons, sons2 in delta_min_sons:
						sp, spp = sons
						sp2, spp2 = sons2
						inds = find_all_mins([c[gp.nr][s.nr]+inn[gpp.nr][spp2.nr]+(numOf_ssons-1)*Ploss,
										c[gp.nr][s.nr]+inn[gpp.nr][sp2.nr]+(numOf_ssons-1)*Ploss,
										c[gpp.nr][s.nr]+inn[gp.nr][spp.nr]+(numOf_ssons-1)*Ploss,
										c[gpp.nr][s.nr]+inn[gp.nr][sp.nr]+(numOf_ssons-1)*Ploss,
										c[gp.nr][s.nr]+c[gpp.nr][s.nr],
										inn[gp.nr][sp.nr]+inn[gpp.nr][spp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss,
										inn[gp.nr][spp.nr]+inn[gpp.nr][sp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss,
										inn[gp.nr][sp.nr]+inn[gpp.nr][sp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss,
										inn[gp.nr][spp.nr]+inn[gpp.nr][spp2.nr]+(numOf_ssons-1)*Ploss+(numOf_ssons-1)*Ploss])

						for ind in inds:
							if ind == 0: 
								if not (-1,spp2.nr) in added_pairs: 
									gppmaplist = unknownleaves_rec[gpp.nr][spp2.nr]  if "?" in gpp.label else inn_rec[gpp.nr][spp2.nr]
									
									for gppmap in gppmaplist:
										reconstruct_delta.append(((gp.nr,s.nr),(gpp.nr,gppmap), (gfath,sfath),"D") )
									added_pairs[(-1, spp2.nr)] = 0
							elif ind == 1: 
								if not (-1,sp2.nr) in added_pairs: 
									gppmaplist = unknownleaves_rec[gpp.nr][sp2.nr]  if "?" in gpp.label else inn_rec[gpp.nr][sp2.nr]
									added_pairs[(-1, sp2.nr)] = 0
									for gppmap in gppmaplist:
										reconstruct_delta.append(((gp.nr,s.nr),(gpp.nr,gppmap), (gfath,sfath),"D") )
							elif ind == 2: 
								if not (spp.nr,-1) in added_pairs: 
									gpmaplist = unknownleaves_rec[gp.nr][spp.nr]  if "?" in gp.label else inn_rec[gp.nr][spp.nr]
									added_pairs[(spp.nr, -1)] = 0
									for gpmap in gpmaplist:
										reconstruct_delta.append(((gpp.nr,s.nr),(gp.nr,gpmap), (gfath,sfath),"D") )
							elif ind == 3: 
								if not (sp.nr,-1) in added_pairs: 
									gpmaplist = unknownleaves_rec[gp.nr][sp.nr]  if "?" in gp.label else inn_rec[gp.nr][sp.nr]
									added_pairs[(sp.nr, -1)] = 0
									for gpmap in gpmaplist:
										reconstruct_delta.append(((gpp.nr,s.nr),(gp.nr,gpmap), (gfath,sfath),"D") )
							elif ind == 4: 
								if not (-1,-1) in added_pairs: 
									added_pairs[(-1, -1)] = 0
									reconstruct_delta.append(((gp.nr,s.nr),(gpp.nr,s.nr), (gfath,sfath),"D"))
							elif ind == 5: 
								if not (sp.nr,spp2.nr) in added_pairs: 
									added_pairs[(spp.nr, spp2.nr)] = 0
									gpmaplist = unknownleaves_rec[gp.nr][sp.nr]  if "?" in gp.label else inn_rec[gp.nr][sp.nr]
									gppmaplist = unknownleaves_rec[gpp.nr][spp2.nr]  if "?" in gpp.label else inn_rec[gpp.nr][spp2.nr]
									for mapspair in list(itertools.product(gpmaplist, gppmaplist)):
										reconstruct_delta.append(((gp.nr,mapspair[0]),(gpp.nr,mapspair[1]), (gfath,sfath),"D"))
							elif ind == 6: 
								if not (spp.nr,sp2.nr) in added_pairs: 
									added_pairs[(spp.nr, sp2.nr)] = 0
									gpmaplist = unknownleaves_rec[gp.nr][spp.nr]  if "?" in gp.label else inn_rec[gp.nr][spp.nr]
									gppmaplist = unknownleaves_rec[gpp.nr][sp2.nr]  if "?" in gpp.label else inn_rec[gpp.nr][sp2.nr]
									for mapspair in list(itertools.product(gpmaplist, gppmaplist)):
										reconstruct_delta.append(((gp.nr,mapspair[0]),(gpp.nr,mapspair[1]), (gfath,sfath),"D"))
							elif ind == 7: 
								if not (sp.nr,sp2.nr) in added_pairs: 
									added_pairs[(sp.nr, sp2.nr)] = 0
									gpmaplist = unknownleaves_rec[gp.nr][sp.nr]  if "?" in gp.label else inn_rec[gp.nr][sp.nr]
									gppmaplist = unknownleaves_rec[gpp.nr][sp2.nr]  if "?" in gpp.label else inn_rec[gpp.nr][sp2.nr]
									for mapspair in list(itertools.product(gpmaplist, gppmaplist)):
										reconstruct_delta.append(((gp.nr,mapspair[0]),(gpp.nr,mapspair[1]), (gfath,sfath),"D"))
							elif ind == 8: 
								if not (spp.nr,spp2.nr) in added_pairs: 
									added_pairs[(spp.nr, spp2.nr)] = 0
									gpmaplist = unknownleaves_rec[gp.nr][spp.nr]  if "?" in gp.label else inn_rec[gp.nr][spp.nr]
									gppmaplist = unknownleaves_rec[gpp.nr][spp2.nr]  if "?" in gpp.label else inn_rec[gpp.nr][spp2.nr]
									for mapspair in list(itertools.product(gpmaplist, gppmaplist)):
										reconstruct_delta.append(((gp.nr,mapspair[0]),(gpp.nr,mapspair[1]), (gfath,sfath),"D"))

					reconstruct_theta = []
					if s.nr != 0: 
						c_theta[g.nr][s.nr] = Ptheta + min(inn[gp.nr][s.nr] + out[gpp.nr][s.nr], inn[gpp.nr][s.nr] + out[gp.nr][s.nr])
						inds = find_all_mins([Ptheta + inn[gp.nr][s.nr] + out[gpp.nr][s.nr],Ptheta + inn[gpp.nr][s.nr] + out[gp.nr][s.nr]])
						for ind in inds:
							if ind == 0: 
								for out_result in out_res[gpp.nr][s.nr]:
									reconstruct_theta.append(((gp.nr,s.nr),(gpp.nr,"->",out_result), (gfath,sfath),"T"))
							else:
								for out_result in out_res[gp.nr][s.nr]:
									reconstruct_theta.append(((gpp.nr,s.nr),(gp.nr,"->", out_result), (gfath,sfath),"T"))

					c[g.nr][s.nr] = min(c_delta[g.nr][s.nr], c_theta[g.nr][s.nr], c_sigma[g.nr][s.nr])
					inn[g.nr][s.nr] = min([c[g.nr][s.nr]] + [inn[g.nr][son.nr] + (numOf_ssons-1)*Ploss for son in s.sons])
					inds=find_all_mins([c[g.nr][s.nr]] + [inn[g.nr][son.nr] + (numOf_ssons-1)*Ploss for son in s.sons])
					for ind in inds:
						if ind == 0:
							inn_rec[g.nr][s.nr].append(s.nr)
						else:
							inn_rec[g.nr][s.nr].extend(inn_rec[g.nr][s.sons[ind-1].nr])
					inAlt[g.nr][s.nr] = min([c[g.nr][s.nr]] + [inAlt[g.nr][son.nr] for son in s.sons])
					inds = find_all_mins([c[g.nr][s.nr]] + [inAlt[g.nr][son.nr] for son in s.sons])
					for ind in inds:
						if ind == 0:
							inAlt_rec[g.nr][s.nr].append(s.nr)
						else:
							inAlt_rec[g.nr][s.nr].extend(inAlt_rec[g.nr][s.sons[ind-1].nr])
					
					inds = find_all_mins([c_delta[g.nr][s.nr], c_theta[g.nr][s.nr], c_sigma[g.nr][s.nr]])
					events[g.nr][s.nr] = ind
					for ind in inds:
						if ind == 0: 
							reconstruct[g.nr][s.nr].extend(reconstruct_delta)
						elif ind == 1: 
							reconstruct[g.nr][s.nr].extend(reconstruct_theta)
						else: 
							reconstruct[g.nr][s.nr].extend(reconstruct_sigma)
			if gp.leaf: cos3(gp, s)
			if gpp.leaf: cos3(gpp, s)	
			cos2(g,s)
			
			cos3(g, s)		
	cos(G, S)
	

	min_inds = find_all_mins(c[G.nr][:-1].tolist())
	for ind in min_inds: count_results(0, ind)
	
	optimal_root_mappings = {sum(k):[] for k in result_counters[0]}
	
	for i in xrange(len(result_counters[0])-1): 
		optimal_root_mappings[sum(result_counters[0][i])].append(i)
	sorted_weights = sorted(optimal_root_mappings.keys())
	sorted_weights.remove(0)
	summ = sum(optimal_root_mappings.keys())

	def random_scenario():
		
		i = random.random() * summ
		maxweight, wsum = 0, 0
		for w in sorted_weights:
			wsum += w
			if i < wsum:
				maxweight = w
				break
		ind=random.randint(0, len(optimal_root_mappings[maxweight])-1)
		chosen_min_rooting_ind = optimal_root_mappings[maxweight][ind]
	
		random_result(0,chosen_min_rooting_ind)

	print "duplication =", Pdelta
	print "loss =", Ploss
	print "HGT =", Ptheta

	print "number of optimal results", sum([sum(k) for k in result_counters[0]])
	print "minimum", min(c[G.nr])
	print 
	
	


	num_of_tests = 1000
	for i in xrange(num_of_tests):
		random_scenario()


	for g in heatmap_results:
		nodeG = GnodesByID[g]
		print nodeG.label+":\t",
		for s in heatmap_results[g]:
			if heatmap_results[g][s]: 
				nodeS = SnodesByID[s]
				print nodeS.label+" "+str(heatmap_results[g][s]/float(num_of_tests))+"\t",
		print

	

