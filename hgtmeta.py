import numpy
import copy
import random
import itertools
import sys
from tree import Node
import getopt

def read(cluster, tree, i):
	while i < len(cluster):
		x = cluster[i]
		if x == "(":
			newtree, i = read(cluster, [], i+1)
			tree.append(newtree)
			i+=1
		elif x == ")":
			return tree, i
		elif x == "," or x == ";" or x == " ":
			i+=1
		else:				
			s = ""
			while x != "," and x != "]" and x != "[" and x != ")" and x != "(" and x != ";":
				s+=x
				i+=1
				x = cluster[i]	
			ss = s.split(":")[0] #removing branch lenghts
			if ss != "": tree.append(ss)
	return tree, i


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


G.lca(S)


leaves_G = G.leaves()
leaves_S = S.leaves()


niesk = float("inf")

c=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
c_sigma=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
c_delta=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
c_theta=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
inn=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
inAlt=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
out=numpy.array([[niesk]*(size_S+1)]*(size_G+1))
out_res = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]
events = numpy.array([['-']*(size_S+1)]*(size_G+1))
reconstruct = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]
result_counters = [[[] for i in  xrange(size_S+1)] for j in xrange(size_G+1)]





Result = [[]]
	
def recontruct(g, s, recon, G, S):
	nodeG = G.nodeById(g)
	nodeS = S.nodeById(s)
	mappings = recon[g][s]
	if len(mappings) == 1:
		mp = mappings[0]
	
		if mp:
			if typee == "M":
				for res in Result:  res.append("%s (%d) -> %s (%d) " %(nodeG.label, g, nodeS.label, s))
				recontruct(mp[0][0], mp[0][1], recon, G, S) 

			else:	
				for res in Result:  res.append("%s (%d) -> %s (%d) transfer to edge (%d, %d)" %(nodeG.label, g, nodeS.label, s, mp[4][0], mp[4][1]))
			if mp[3] == "D" or mp[3] == "S":
				for res in Result:  res.append("%s (%d) -> %s (%d) %s" %(nodeG.label, g, nodeS.label, s, mp[3]))
				recontruct(mp[0][0], mp[0][1], recon, G, S) 
				recontruct(mp[1][0], mp[1][1], recon, G, S)
			else:
				recontruct(mp[0][0], mp[0][1], recon, G, S)
				recontruct(mp[1][0], mp[1][1], recon, G, S)  
					
			
		else:
			for res in Result:  res.append("%s (%d) -> %s (%d)" %(nodeG.label, g, nodeS.label, s))

	else:
		for res in Result:  res.append("%s (%d) -> %s (%d)" %(nodeG.label, g, nodeS.label, s))


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
				numOfCh = sum(result_counters[res[0][0]][res[0][-1]])+ sum(result_counters[res[1][0]][res[1][-1]])
				if numOfCh == 0: 
					result_counters[X][Y].append(1)
				else: 
					result_counters[X][Y].append(numOfCh)

			
		else: 
			result_counters[X][Y].append(0)

heatmap_results = {k.nr:{kk.nr:0 for kk in S.leaves()+S.inner()} for k in G.leaves()}	
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
	nodeG = G.nodeById(g)
	nodeS = S.nodeById(s)
	if not nodeS.leaf:
		if not "?" in nodeG.label:
			nodeS = nodeG.map
		else:
			while not ss.leaf:
				i = random.randint(len(ss.sons)-1)
				ss = ss.sons[i]
			nodeS = ss

	
	if nodeG.leaf: 
		heatmap_results[nodeG.nr][nodeS.nr] += 1
	if len(result) >= 1:
		random_result(res[0][0], res[0][-1])
		random_result(res[1][0], res[1][-1])
	

if True:	
	for g in leaves_G:
		if "?" in g.label:
			for i in xrange(len(c[g.nr])-1): 
				#print i
				nodeS = S.nodeById(i)
				if nodeS.leaf: c[g.nr][i] = 0
		else: c[g.nr][g.map.nr] = 0

		if "?" not in g.label:
			inn[g.nr][g.map.nr] = Ploss*0 
			inAlt[g.nr][g.map.nr] = 0
			ps = g.map.fath
			dS = 1#path between s and Le(g)
			while ps != None:
				inn[g.nr][ps.nr] = Ploss*dS
				inAlt[g.nr][ps.nr] = 0
	    			ps=ps.fath
				dS += 1
		else:
			for gmap in S.leaves():
				ps = gmap.fath
				inn[g.nr][gmap.nr] = Ploss*0
				inAlt[g.nr][gmap.nr] = 0
				dS = 1
				while ps !=None:
					inn[g.nr][ps.nr] = Ploss*dS
					inAlt[g.nr][ps.nr] = 0
					ps = ps.fath
					dS+=1


	g=G
	Fi = {}
	def cos(g, s):
		if not g.leaf:
			gp, gpp = g.sons
			cos(gp, s) 
			cos(gpp, s) 	
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
									s_siblings = [son for son in s.fath.sons if son.nr != s.nr]
									for s_sib in s_siblings:
										if not s_sib.leaf:
											minn = min([inn[g.nr][s_sib_son.nr] for s_sib_son in s_sib.sons])
											inds2 = find_all_mins([inn[g.nr][s_sib_son.nr] for s_sib_son in s_sib.sons])
											for ind2 in inds2:
												rec.append(s_sib.sons[ind2].nr)
										else: rec.append(s_sib.nr)
								else:
									rec.append(sson.nr)
							else:
								spson = other_ssons[ind-1]
								if not spson.leaf:
									minn = min([inn[g.nr][spson_son.nr] for spson_son in spson.sons])	
									inds2 = find_all_mins([inn[g.nr][spson_son.nr] for spson_son in spson.sons])	
									for ind2 in inds2:
										rec.append(spson.sons[ind2].nr)	
								else: rec.append(spson.nr)
						out[g.nr][sson.nr] = out_g_sson
						out_res[g.nr][sson.nr] = copy.deepcopy(rec)


					for sson in s.sons: cos3(g, sson)
									
			def cos2(g, s):		
				gfath = g.fath.nr if g.nr != 0 else 0
				sfath = g.fath.map.nr if g.nr != 0 else 0
				if s.leaf:
					c_sigma[g.nr][s.nr] = niesk
					c_delta[g.nr][s.nr] = Pdelta + c[gp.nr][s.nr] + c[gpp.nr][s.nr]
					reconstruct_delta = [((gp.nr,s.nr),(gpp.nr,s.nr), (gfath,sfath), "D")]
					if s.nr != 0: #if it's not root
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
					inAlt[g.nr][s.nr] = c[g.nr][s.nr]

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
					numOf_sigma_losses = numOf_ssons - 2 # -2 bc g is binary and it is number of losses
					for sp, spp in ssons_pairs:
						minn = min(inn[gp.nr][sp.nr] + inn[gpp.nr][spp.nr], inn[gpp.nr][sp.nr] + inn[gp.nr][spp.nr])
						if minn < sigma_min: 
							sigma_min = minn
							sigma_min_sons = [(sp, spp)]
						elif minn == sigma_min: sigma_min_sons.append((sp, spp))
					reconstruct_sigma = []
					c_sigma[g.nr][s.nr] = sigma_min + Ploss * numOf_sigma_losses
					for sp, spp in sigma_min_sons:
						
						inds = find_all_mins([inn[gp.nr][sp.nr] + inn[gpp.nr][spp.nr], inn[gpp.nr][sp.nr] + inn[gp.nr][spp.nr]])
						for ind in inds:
							if ind == 0: reconstruct_sigma.append(((gp.nr,sp.nr),(gpp.nr,spp.nr), (gfath,sfath),"S"))
							else: reconstruct_sigma.append(((gpp.nr,sp.nr),(gp.nr,spp.nr),(gfath,sfath),"S"))
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
									reconstruct_delta.append(((gp.nr,s.nr),(gpp.nr,spp2.nr), (gfath,sfath),"D") )
									added_pairs[(-1, spp2.nr)] = 0
							elif ind == 1: 
								if not (-1,sp2.nr) in added_pairs: 
									added_pairs[(-1, sp2.nr)] = 0
									reconstruct_delta.append(((gp.nr,s.nr),(gpp.nr,sp2.nr),(gfath,sfath),"D") )
							elif ind == 2: 
								if not (spp.nr,-1) in added_pairs: 
									added_pairs[(spp.nr, -1)] = 0
									reconstruct_delta.append(((gpp.nr,s.nr),(gp.nr,spp.nr), (gfath,sfath),"D") )
							elif ind == 3: 
								if not (sp.nr,-1) in added_pairs: 
									added_pairs[(sp.nr, -1)] = 0
									reconstruct_delta.append(((gpp.nr,s.nr),(gp.nr,sp.nr), (gfath,sfath),"D"))
							elif ind == 4: 
								if not (-1,-1) in added_pairs: 
									added_pairs[(-1, -1)] = 0
									reconstruct_delta.append(((gp.nr,s.nr),(gpp.nr,s.nr), (gfath,sfath),"D"))
							elif ind == 5: 
								if not (sp.nr,spp2.nr) in added_pairs: 
									added_pairs[(spp.nr, spp2.nr)] = 0
									reconstruct_delta.append(((gp.nr,sp.nr),(gpp.nr,spp2.nr), (gfath,sfath),"D"))
							elif ind == 6: 
								if not (spp.nr,sp2.nr) in added_pairs: 
									added_pairs[(spp.nr, sp2.nr)] = 0
									reconstruct_delta.append(((gp.nr,spp.nr),(gpp.nr,sp2.nr), (gfath,sfath),"D"))
							elif ind == 7: 
								if not (sp.nr,sp2.nr) in added_pairs: 
									added_pairs[(sp.nr, sp2.nr)] = 0
									reconstruct_delta.append(((gp.nr,sp.nr),(gpp.nr,sp2.nr), (gfath,sfath),"D"))
							elif ind == 8: 
								if not (spp.nr,spp2.nr) in added_pairs: 
									added_pairs[(spp.nr, spp2.nr)] = 0
									reconstruct_delta.append(((gp.nr,spp.nr),(gpp.nr,spp2.nr), (gfath,sfath),"D"))


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
					inAlt[g.nr][s.nr] = min([c[g.nr][s.nr]] + [inAlt[g.nr][son.nr] for son in s.sons])

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
	print

	for g in heatmap_results:
		nodeG = G.nodeById(g)
		print nodeG.label+":\t",
		for s in heatmap_results[g]:
			if heatmap_results[g][s]: 
				nodeS = S.nodeById(s)
				print nodeS.label+" "+str(heatmap_results[g][s]/float(num_of_tests))+"\t",
		print

