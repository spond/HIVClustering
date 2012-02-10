import csv, argparse, sys, datetime, time, random, os.path
from scipy import stats
from BeastNetwork import *
import _hyphyinterface as hy

#-------------------------------------------------------------------------------		

def import_edi (file):
	edi_by_id = {}
	ediReader = csv.reader(file)
	header = ediReader.next ()
	if len (header) != 14:
		raise Exception ('Expected a .csv file with 14 columns as input')
	
	for line in ediReader:
		if len (line[1]): # has PID
			id = line[1].replace ('-','')
		else:
			id = line[0]
		
		geno_date = None
		if len (line[2]): # geno
			geno_date = time.strptime (line[2],'%m/%d/%Y')
		
		
		drug_date = None
		if len (line[4]): # drugz
			drug_date = time.strptime (line[4],'%m/%d/%Y')
			
		edi_date = None
		if len (line [6]): # edi
			edi_date = time.strptime (line[6],'%m/%d/%Y')
			
		naive 	 = False
		if line[3] == 'ARV Naive':
			naive = True
			
		if geno_date and edi_date:
			if edi_date > geno_date:
				#print time.mktime(edi_date) - time.mktime(geno_date)
				
				part1 = time.strftime ("%m/%d",edi_date)
				part2 = time.strftime ("%Y",geno_date)
				new_edi_date = time.strptime ("/".join((part1,part2)),'%m/%d/%Y')
				#edi_date.tm_year = geno_date.tm_year
				if new_edi_date > geno_date:
					continue
				else:	
					edi_date = new_edi_date
				
		viral_load = None
		if len (line [8]): # vl
			viral_load = int (line[8])
			
		edi_by_id [id] = [geno_date, drug_date, edi_date, viral_load, naive]
		#if (edi_date and drug_date and edi_date > drug_date):
		#	print "Fail %s" % id, edi_date, drug_date
		
	return edi_by_id
	
#-------------------------------------------------------------------------------		

def split_nodes_by_year (nodes):
	by_year = {}
	for node in nodes:
		year = node.get_baseline_date ()
		if year not in by_year:
			by_year [year] = [node]
		else:
			by_year [year].append (node)
	return by_year

#-------------------------------------------------------------------------------		

def match_sample_by_year (case, control):
	by_year_case    = split_nodes_by_year (case)
	by_year_control = split_nodes_by_year (control)
	
	matched_sample  = []
	
	for k in by_year_case:
		for k2 in random.sample(by_year_control[k], min(len(by_year_case[k]),len(by_year_control[k]))):
			matched_sample.append (k2)
			
	return matched_sample
	
#-------------------------------------------------------------------------------		

def degree_to_list (deg,shift = False):
	listed = []
	for k in range (len(deg)):
		for z in range (deg[k]):
			listed.append(k + (1 if shift else 0))
	return listed
	
#-------------------------------------------------------------------------------		

def compare_degrees (deg1, deg2, do_convert = True, shift = False):
	if (do_convert):
		d1 = degree_to_list(deg1, shift)
		d2 = degree_to_list(deg2, shift)
	else:
		d1 = deg1
		d2 = deg2
	
	print d1
	print d2
	
	print sum(d1)/float(len(d1)), sum(d2)/float(len(d2))
	return stats.mannwhitneyu (d1,d2)
	

#-------------------------------------------------------------------------------		

random.seed ()

arguments = argparse.ArgumentParser(description='Read filenames.')

arguments.add_argument('-i', '--input',   help = 'Input CSV file with inferred genetic links (or stdin if omitted). Must be a CSV file with three columns: ID1,ID2,distance.')
arguments.add_argument('-d', '--dot',   help = 'Output DOT file for GraphViz (or stdout if omitted)')
arguments.add_argument('-c', '--cluster', help = 'Output a CSV file with cluster assignments for each sequence', required = True)
arguments.add_argument('-t', '--threshold', help = 'Only count edges where the distance is less than this threshold')
arguments.add_argument('-e', '--edi',   help = 'A .csv file with EDI dates')
arguments.add_argument('-f', '--format',   help = 'Sequence ID format. One of AEH (ID | sample_date | otherfiels default), LANL (e.g. B_HXB2_K03455_1983 : subtype_country_id_year -- could have more fields), plain (treat as sequence ID only, no meta)')
arguments.add_argument('-x', '--exclude',   help = 'Exclude any sequence which belongs to a cluster containing a "reference" strain, defined by the year of isolation. The value of this argument is an integer year (e.g. 1983) so that any sequence isolated in or before that year (e.g. <=1983) is considered to be a lab strain. This option makes sense for LANL or AEH data.')

settings = arguments.parse_args()

if settings.input == None:
	settings.input = sys.stdin
else:
	try:
		settings.input = open (settings.input, 'r')
	except IOError:
		print "Failed to open '%s' for reading" % (settings.input)
		raise
	
if settings.dot == None:
	settings.dot = sys.stdout
else:
	try:
		settings.dot = open (settings.dot, 'w')
	except IOError:
		print "Failed to open '%s' for writing" % (settings.dot)
		raise
		
	
edi = None

if settings.edi != None:
	try:
		settings.edi = open (settings.edi, 'r')
		edi = import_edi (settings.edi)
	except IOError:
		print "Failed to open '%s' for reading" % (settings.edi)
		raise

try:
	settings.cluster = open (settings.cluster, 'w')
except IOError:
	print "Failed to open '%s' for writing" % (settings.cluster)
	raise
	
formatter = parseAEH

if settings.format != None:
	formats = {"AEH" : parseAEH, "LANL": parseLANL, "plain": parsePlain}
	try:
		formatter = formats[settings.format]
	except KeyError:
		print "%s is not a valid setting for 'format' (must be in %s)" % (settings.edi,str(formats.keys()))
		raise
		
if settings.exclude != None:
	try:
		settings.exclude = datetime.datetime (int (settings.exclude), 12, 31)	
	except ValueError:
		print "Invalid contaminant threshold year '%s'" % (settings.exclude)
		raise 

edgeReader = csv.reader(settings.input)

if settings.threshold == None:
	settings.threshold = 1000.0
else:
	settings.threshold = float (settings.threshold)

header = edgeReader.next ()

if len (header) != 3:
	raise Exception ('Expected a .csv file with 3 columns as input')
	
nodesRead = {}	
degrees   = {}


network = transmission_network ()

for line in edgeReader:
	distance = float(line[2])
	
	if distance > settings.threshold:
		continue

	network.add_an_edge(line[0],line[1],distance,formatter)	
	 
network.compute_clusters()
clusters = network.retrieve_clusters ()

if settings.exclude != None:
	remove_these_ids = []
	for c in clusters:
		for node in clusters[c]:
			if tm_to_datetime (node.dates[0]) <= settings.exclude:
				if node.cluster_id not in remove_these_ids:
					remove_these_ids.append(node.cluster_id)
		
	print "Removing %d clusters with contaminant sequences" % len (remove_these_ids)
	network.apply_cluster_filter (remove_these_ids)
	network.compute_clusters()
	clusters = network.retrieve_clusters ()
	
network_stats = network.get_edge_node_count ()
print "Read a graph on %d nodes with %d edges" % (network_stats['nodes'], network_stats['edges'])
print "Found %d clusters" % len(clusters)
print "Maximum cluster size = %d nodes" % max ([len (clusters[c]) for c in clusters])


network.write_clusters (settings.cluster)

print "Fitting the degree distribution to various densities"
degrees = network.get_degree_distribution()



if settings.dot != None:
	network.generate_dot (settings.dot)

hy_instance = hy.HyphyInterface ();
	
script_path = os.path.realpath(__file__)
hbl_path =  os.path.join(os.path.dirname(script_path), "HBL", "DegreeDistributions.bf")

all_deg = network.get_degree_distribution()
hy_instance.queuevar ('allDegs', all_deg)
hy_instance.runqueue (batchfile = hbl_path)

#print hy_instance.stdout
bestDistro = hy_instance.getvar ('BestDistro',hy.HyphyInterface.STRING)
bestRho	   = hy_instance.getvar (bestDistro, hy.HyphyInterface.NUMBER);

print "Best distribution is '%s' with rho = %g" % (bestDistro, bestRho)

