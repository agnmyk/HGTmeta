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

class Node:
	w=0
	def __init__(self, n, f=None, d=0):
		self.nr = Node.w
		Node.w+=1
		self.fath = f
		self.map = None
		self.deep = d
		self.time = 0
		if type(n)==str:
			self.leaf = True
			self.label = n
			
		else:
			self.leaf = False
			self.sons = [Node(s, self, d+1) for s in n]
			self.label = "n%d" %self.nr

	@classmethod
	def readTree(cls,src):
		q=read(src.rstrip(), [], 0)[0][0]
		cls.w=0
		S=cls(q)
		return S 

	def leaves(self):
		if self.leaf:	return [self]
		else: 
			l = []
			for s in self.sons: l+=s.leaves()
			return l

	def inner(self):
		if not self.leaf: 
			l = [self]
			for s in self.sons: 
				if not s.leaf: l+=s.inner()
			return l

	def lca(self, node):
		if self.leaf:
			for i in node.leaves():		
				if i.label==self.label or "?" in self.label:
					self.map=i
					break
		else:
			for s in self.sons: 
				if s.map==None: s.lca(node)
			r = [s.map for s in self.sons]
			
			min_r = min([(rr, rr.deep) for rr in r], key=lambda x:x[1])[0]
			for rr in r:
				while rr.deep > min_r.deep:
					rr = rr.fath
				while rr != min_r:
					rr = rr.fath
					min_r= min_r.fath
			self.map=min_r
			


	def printTree(self):
		if self.leaf:
			return self.label + "#" + str(self.nr)
		else:
			return '('+' '.join([s.printTree() for s in self.sons])+')'

	def printMap(self):
		print str(self.label)+'->'+str(self.map.label)
		if not self.leaf:
			for s in self.sons: s.printMap()



	def printTime(self):
		if not self.leaf: 
			self.sons[0].printTime()
			self.sons[1].printTime()
		print "%d, %d" %(self.nr, self.time)

	def nodeById(self, nr):
		if self.nr == nr: return self
		if not self.leaf:
			node = None
			c = 0
			while not node and c < len(self.sons):
				node = self.sons[c].nodeById(nr)
				c += 1
		else: return None
		return node
