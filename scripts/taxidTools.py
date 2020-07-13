#!/usr/bin/env python3
"""
Provides the Taxdump class to load ncbi taxonomy dumpfiles and work with taxids and lineages. 
"""


__version__ = "1.2.1"
__author__ = "Gregoire Denay"
__email__ = "gregoire.denay@cvua-rrw.de"


from collections import UserDict

FULL_TAXONOMY = ('forma',
				 'varietas',
				 'subspecies',
				 'species',
				 'species subgroup', 
				 'species group', 
				 #'series', Inconsistent use between botany and zoology
				 #'subsection', Inconsistent use between botany and Zoology
				 #'section', Inconsistent use between botany and Zoology
				 'subgenus',
				 'genus',
				 'subtribe',
				 'tribe',
				 'subfamily',
				 'family',
				 'superfamily',
				 'parvorder',
				 'infraorder',
				 'suborder',
				 'order',
				 'superorder',
				 'subcohort',
				 'cohort',
				 'infraclass',
				 'subclass'
				 'class',
				 'superclass',
				 'subphylum',
				 'phylum',
				 'superphylum',
				 'subkingdom',
				 'kingdom',
				 'superkingdom',
				 )
				 
LINNE_TAXONOMY = ('species', 
				  'genus', 
				  'family', 
				  'order', 
				  'class', 
				  'phylum', 
				  'kingdom', 
				  'superkingdom',
				  )


class ReadOnlyDict(UserDict):
	"""
	Readonly subclass of UserDict. 
	Can still be updated by using ReadOnlyDict.data[key] = val (Should only be used for object initialization)
	"""
	def __init__(self, mapping=None):
		"""
		A subclass of Userdict without self-modification capabilities.
		
		Arguments:
		----------
		mapping: dict
			Dict to populate the object
		"""
		super().__init__()
		if mapping:
			for key, val in dict(mapping).items():
				self.data[key] = val
		
	def __readonly__(self, *args, **kwargs):
		"raise NotImplementedError"
		raise NotImplementedError("This object is read-only")
		
	__setitem__ = __readonly__
	__delitem__ = __readonly__
	pop = __readonly__
	popitem = __readonly__
	clear = __readonly__
	update = __readonly__
	setdefault = __readonly__
	del __readonly__
	
	
class Taxdump(object):
	"""
	Reads in taxdump files and provides methods to work on lineages.
	
	Methods returning taxonomic lineages will only return ranks based on the currently in use taxonomy. This is either defined by passing a list of
	ranks during object creation, modifying it through the setTaxonomy() method or directly as an argument of the methods.
	"""
	def __init__(self, rankedlineage_dmp, nodes_dmp, want_ranks=None):
		"""
		Create a Taxdump object.
		
		Arguments
		---------
		rankedlineage_dmp: str
			Path to rankedlineage.dmp
		nodes_dmp: str
			Path to nodes.dmp
		want_ranks: list or tuple 
			Iterator of ranks as defined by the NCBI taxonomy nomenclature. Ordered from the lowest to the highest rank.
		"""
		# useful data will be stored in a dict as {taxid: (name, rank, parent_taxid)} 
		# the aim here is to keep data storage at a minimum to use as little memory as nescessary without the hurdle of indexing the files
		self._data = ReadOnlyDict()
		
		# Parsing name file
		tmp_name = {} # temporary storage for taxid-name pairs
		for line in _parse_dump(rankedlineage_dmp):
			# retriving unique name for each taxid
			taxid = line[0]
			name = line[1]
			tmp_name[taxid] = name

		# Parsing node file and creating the _data dict on the fly
		for line in _parse_dump(nodes_dmp):
			taxid = line[0]
			name = tmp_name.pop(taxid)
			rank = line[2]
			parent = line[1]
			self._data.data[taxid] = (name, rank, parent)
			
		# defining taxonomy 
		self.taxonomy_full = FULL_TAXONOMY
		self.taxonomy_linne = LINNE_TAXONOMY
		
		if want_ranks:
			self.want_ranks = want_ranks
		else:
			self.want_ranks = self.taxonomy_linne
	
	def taxonomyInUse(self):
		"Display the taxonomic ranks currently in use"
		return self._taxonomyView(self.want_ranks)
	
	def setTaxonomy(self, ranks):
		"""
		Sets the default taxonomy ranks
		
		ranks: list or tuple
			Iterator of ranks from the lowest rank to the highest
		"""
		self.want_ranks= ranks
		
	def getName(self, taxid):
		"""
		Get taxid name
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
			
		Returns
		-------
		str:
			node name
		"""
		return self._data[taxid][0]
		
	def getRank(self, taxid):
		"""
		Get taxid rank
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
			
		Returns
		-------
		str:
			node rank
		"""
		return self._data[taxid][1]		
	
	def getParent(self, taxid):
		"""
		Retrieve parent taxid 
		
		Arguments:
		----------
		taxid: str
			Taxonomic identification number
		
		Returns
		-------
		str:
			Parent taxid
		"""
		return self._data[taxid][2]
		
	def getLineageAsDict(self, taxid, want_ranks=None, asNames=False):
		"""
		Get taxid ancestry
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
		want_ranks: list
			List of wanted taxonomic ranks from the lowest node to the highest node
		asNames: bool
			Report taxonomy with node names (default reports taxid)
		
		Returns
		-------
		dict:
			Dictionnary of taxonomy ranks and values.
			Note that missing ranks will be omitted.
		"""
		want_ranks = want_ranks or self.want_ranks
		
		fulllineage = self.getFullLineage(taxid)
		d = self._filterLineage(fulllineage, want_ranks)
		
		if asNames:
			for key, val in d.items():
				d[key]=self.getName(val)
				
		return d
		
	def getLineageAsList(self, taxid, want_ranks=None, asNames=False):
		"""
		Get taxid ancestry
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
		want_ranks: list
			List of wanted taxonomic ranks from the lowest node to the highest node
		asNames: bool
			Report taxonomy with node names (default reports taxid)
		
		Returns
		-------
		list:
			List of taxid values from the lowest to the highest node.
			Note that missing ranks will be returned as an empty string
		"""
		want_ranks = want_ranks or self.want_ranks
		
		l = []
		lineage = self.getLineageAsDict(taxid, want_ranks=want_ranks, asNames=asNames)
		for rank in want_ranks:
			l.append(lineage[rank])
		return l
		
	def getLineageAsView(self, taxid, want_ranks=None):
		"""
		Display taxid ancestry
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
		want_ranks: list
			List of wanted taxonomic ranks from the lowest node to the highest node
		asNames: bool
			Report taxonomy with node names (default reports taxid)
		"""
		want_ranks = want_ranks or self.want_ranks
		
		l=[]
		for key, val in self.getLineageAsDict(taxid, want_ranks=want_ranks).items():
			l.append("{} ({}: {})".format(self.getName(val), key, val))
		return self._taxonomyView(l)
		
	def lowestCommonNode(self, taxid_list, want_ranks=None):
		"""
		Get lowest common node of a bunch of taxids
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
		want_ranks: list
			List of wanted taxonomic ranks from the lowest node to the highest node
			
		Returns
		-------
		str: 
			Lowest common taxid
		"""
		want_ranks = want_ranks or self.want_ranks
		
		# Collect input ranks
		ranks = []
		for taxid in taxid_list:
			rank = self.getRank(taxid)
			if rank in want_ranks:
				ranks.append(rank)
			# If the rank is not in the want rank list, find all the ranks in the lineage and take the lowest one
			else:
				# Getting the list of all ranks in the lineage
				rank_lineage = self.getFullLineage(taxid).keys()
				indices = []
				# finding rank index in want taxonomy
				for r in rank_lineage:
					try:
						indices.append(want_ranks.index(r))
					except ValueError:
						pass
				# recovering lowest rank
				ranks.append(want_ranks[min(indices)])
		
		# Recover lowest possible rank determination
		indices = []
		for rank in set(ranks):
			indices.append(want_ranks.index(rank))
		tax_depth= want_ranks[max(indices):]
		
		# Creating a list of lineage lists, ranks are in the same order
		lineages = []
		for taxid in taxid_list:
			lineages.append(self.getLineageAsList(taxid, tax_depth))
		
		# Use zip to create list of taxid by rank and set to check the number of different entries
		# The rank variable keeps track of the rank index
		rank = 0
		for ranked in zip(*lineages):
			if len(set(tuple(ranked))) == 1: 
				break
			else:
				rank += 1
				
		return lineages[0][rank]
		
	def getFullLineage(self, taxid):
		"""
		Get the full lineage of a taxid
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
		
		Returns
		-------
		dict: 
			Full taxonomic lineage with ranks as keys
		"""
		lineage = {self.getRank(taxid): taxid}
		
		parent = self.getParent(taxid)
		prank = self.getRank(parent)
		
		while parent != '1': # taxid 1 is root 
			lineage[prank] = parent
			parent = self.getParent(parent)
			prank = self.getRank(parent)
		
		return lineage
		
	def getChildren(self, taxid, want_ranks=None, asNames=False):
		"""
		Find the descendants of a taxonomic node
		
		Arguments
		---------
		taxid: str
			Taxonomic identification number
		want_ranks: list
			List of wanted taxonomic ranks from the lowest node to the highest node
		asNames: bool
			Report taxonomy with node names (default reports taxid)
			
		Returns
		-------
		list of dicts:
			List of lineages as a dictionnary. Each downstream node will return a full lineage.
		"""
		want_ranks = want_ranks or self.want_ranks
		
		children = []

		for key in self._data.keys():
			if taxid in self.getFullLineage(key).values():
				lineage = self.getLineageAsDict(key, want_ranks=want_ranks, asNames=asNames)
				children.append(lineage)
				
		return children
		
	def groupByRank(self, taxid_list, rank, silent=True):
		"""
		Cluster a list of taxid at the given rank.
		
		Arguments
		---------
		taxid_list: list of str
			List of taxid to cluster
		rank: str
			rank to which the taxid should be clustered
		silent: bool
			if False, will raise a KeyError if a taxid does not have a parent at the wanted rank.
			if True, taxids that cannot be stratified will be grouped together
			
		Returns
		-------
		dict:
			dict of clusters where keys are the nodes at the clustering rank and values are list of taxid under that node			
		"""
		d= {}
		for taxid in taxid_list:
			lineage = self.getFullLineage(taxid)
			
			# checking if rank is in the lineage
			try: 
				cluster = lineage[rank]
			except KeyError:
				if silent: 
					cluster = "unassigned"
				else: 
					raise KeyError("Taxid: '{}' does not have a parent with the rank '{}'".format(taxid, rank))
			
			# appending taxid to cluster
			try:
				d[cluster].append(taxid)
			except KeyError:
				d[cluster]= [taxid]
		
		return d
		
	@staticmethod	
	def _filterLineage(dict_tax, want_ranks):
		"""
		Filter a dict to return only entries at the wanted ranks. Missing ranks are replaced bvy empty strings.
		"""
		d = {}
		for rank in want_ranks:
			try: 
				d[rank] = dict_tax[rank]
			except KeyError:
				d[rank] = ""
		return d

	@staticmethod
	def _taxonomyView(list_of_ranks):
		"""
		Pretty view of taxonomy
		"""
		print(list_of_ranks[-1])
		if len(list_of_ranks) >1: 
			indent = 1
			for r in list_of_ranks[::-1][1:]:
				print("{}|_{}".format(" "*indent, r))
				indent += 1

	def __len__(self):
		return len(self._data)
		
	def __getitem__(self, index):
		return self._data[index]
		
	def __iter__(self):
		yield from self._data.keys()
	
class GbToTaxid(ReadOnlyDict):
	"""
	Read-only dictionnary to translate Genbank accession numbers to Taxids.
	Can (should?) be used with :
	with GbToTaxid("acc2taxid_file.tsv") as gb2tax:
		# do something
	"""
	def __init__(self, taxidfile):
		"""
		Create a translator
		
		Arguments
		---------
		taxidfile: str
			Path to translation file, assumes that accessions are in the first field and taxid in the second
		"""
		super().__init__(self)
		#print("Parsing {}, be patient...".format(taxidfile))
		with open(taxidfile, 'r') as trans:
			for line in trans:
				l = [item.strip() for item in line.split()]
				self.data[l[0]] = l[1]
		
	def __enter__(self):
		return self
		
	def __exit__(self, exc_type, exc_value, exc_traceback):
		pass
	

def _parse_dump(filepath):
	"""
	Dump file line iterator, returns a list of fields
	"""
	with open(filepath, 'r') as dmp:
		for line in dmp:
			yield [item.strip() for item in line.split("|")]