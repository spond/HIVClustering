import csv, argparse, sys, os.path
import _hyphyinterface as hy
import cStringIO


country_codes = {"AD" : "Andorra", "AE" : "United Arab Emirates", "AF" : "Afghanistan", "AG" : "Antigua and Barbuda", "AI" : "Anguilla", "AL" : "Albania", "AM" : "Armenia", "AO" : "Angola", "AQ" : "Antarctica", "AR" : "Argentina", "AS" : "American Samoa", "AT" : "Austria", "AU" : "Australia", "AW" : "Aruba", "AX" : "land Islands", "AZ" : "Azerbaijan", "BA" : "Bosnia and Herzegovina", "BB" : "Barbados", "BD" : "Bangladesh", "BE" : "Belgium", "BF" : "Burkina Faso", "BG" : "Bulgaria", "BH" : "Bahrain", "BI" : "Burundi", "BJ" : "Benin", "BL" : "Saint Barthlemy", "BM" : "Bermuda", "BN" : "Brunei Darussalam", "BO" : "Bolivia, Plurinational State of", "BQ" : "Bonaire, Sint Eustatius and Saba", "BR" : "Brazil", "BS" : "Bahamas", "BT" : "Bhutan", "BV" : "Bouvet Island", "BW" : "Botswana", "BY" : "Belarus", "BZ" : "Belize", "CA" : "Canada", "CC" : "Cocos (Keeling) Islands", "CD" : "Congo, the Democratic Republic of the", "CF" : "Central African Republic", "CG" : "Congo", "CH" : "Switzerland", "CI" : "Cte d'Ivoire", "CK" : "Cook Islands", "CL" : "Chile", "CM" : "Cameroon", "CN" : "China", "CO" : "Colombia", "CR" : "Costa Rica", "CU" : "Cuba", "CV" : "Cape Verde", "CW" : "Curaao", "CX" : "Christmas Island", "CY" : "Cyprus", "CZ" : "Czech Republic", "DE" : "Germany", "DJ" : "Djibouti", "DK" : "Denmark", "DM" : "Dominica", "DO" : "Dominican Republic", "DZ" : "Algeria", "EC" : "Ecuador", "EE" : "Estonia", "EG" : "Egypt", "EH" : "Western Sahara", "ER" : "Eritrea", "ES" : "Spain", "ET" : "Ethiopia", "FI" : "Finland", "FJ" : "Fiji", "FK" : "Falkland Islands (Malvinas)", "FM" : "Micronesia, Federated States of", "FO" : "Faroe Islands", "FR" : "France", "GA" : "Gabon", "GB" : "United Kingdom", "GD" : "Grenada", "GE" : "Georgia", "GF" : "French Guiana", "GG" : "Guernsey", "GH" : "Ghana", "GI" : "Gibraltar", "GL" : "Greenland", "GM" : "Gambia", "GN" : "Guinea", "GP" : "Guadeloupe", "GQ" : "Equatorial Guinea", "GR" : "Greece", "GS" : "South Georgia and the South Sandwich Islands", "GT" : "Guatemala", "GU" : "Guam", "GW" : "Guinea-Bissau", "GY" : "Guyana", "HK" : "Hong Kong", "HM" : "Heard Island and McDonald Islands", "HN" : "Honduras", "HR" : "Croatia", "HT" : "Haiti", "HU" : "Hungary", "ID" : "Indonesia", "IE" : "Ireland", "IL" : "Israel", "IM" : "Isle of Man", "IN" : "India", "IO" : "British Indian Ocean Territory", "IQ" : "Iraq", "IR" : "Iran, Islamic Republic of", "IS" : "Iceland", "IT" : "Italy", "JE" : "Jersey", "JM" : "Jamaica", "JO" : "Jordan", "JP" : "Japan", "KE" : "Kenya", "KG" : "Kyrgyzstan", "KH" : "Cambodia", "KI" : "Kiribati", "KM" : "Comoros", "KN" : "Saint Kitts and Nevis", "KP" : "Korea, Democratic People's Republic of", "KR" : "Korea, Republic of", "KW" : "Kuwait", "KY" : "Cayman Islands", "KZ" : "Kazakhstan", "LA" : "Lao People's Democratic Republic", "LB" : "Lebanon", "LC" : "Saint Lucia", "LI" : "Liechtenstein", "LK" : "Sri Lanka", "LR" : "Liberia", "LS" : "Lesotho", "LT" : "Lithuania", "LU" : "Luxembourg", "LV" : "Latvia", "LY" : "Libya", "MA" : "Morocco", "MC" : "Monaco", "MD" : "Moldova, Republic of", "ME" : "Montenegro", "MF" : "Saint Martin (French part)", "MG" : "Madagascar", "MH" : "Marshall Islands", "MK" : "Macedonia, the former Yugoslav Republic of", "ML" : "Mali", "MM" : "Myanmar", "MN" : "Mongolia", "MO" : "Macao", "MP" : "Northern Mariana Islands", "MQ" : "Martinique", "MR" : "Mauritania", "MS" : "Montserrat", "MT" : "Malta", "MU" : "Mauritius", "MV" : "Maldives", "MW" : "Malawi", "MX" : "Mexico", "MY" : "Malaysia", "MZ" : "Mozambique", "NA" : "Namibia", "NC" : "New Caledonia", "NE" : "Niger", "NF" : "Norfolk Island", "NG" : "Nigeria", "NI" : "Nicaragua", "NL" : "Netherlands", "NO" : "Norway", "NP" : "Nepal", "NR" : "Nauru", "NU" : "Niue", "NZ" : "New Zealand", "OM" : "Oman", "PA" : "Panama", "PE" : "Peru", "PF" : "French Polynesia", "PG" : "Papua New Guinea", "PH" : "Philippines", "PK" : "Pakistan", "PL" : "Poland", "PM" : "Saint Pierre and Miquelon", "PN" : "Pitcairn", "PR" : "Puerto Rico", "PS" : "Palestinian Territory, Occupied", "PT" : "Portugal", "PW" : "Palau", "PY" : "Paraguay", "QA" : "Qatar", "RE" : "Runion", "RO" : "Romania", "RS" : "Serbia", "RU" : "Russian Federation", "RW" : "Rwanda", "SA" : "Saudi Arabia", "SB" : "Solomon Islands", "SC" : "Seychelles", "SD" : "Sudan", "SE" : "Sweden", "SG" : "Singapore", "SH" : "Saint Helena, Ascension and Tristan da Cunha", "SI" : "Slovenia", "SJ" : "Svalbard and Jan Mayen", "SK" : "Slovakia", "SL" : "Sierra Leone", "SM" : "San Marino", "SN" : "Senegal", "SO" : "Somalia", "SR" : "Suriname", "SS" : "South Sudan", "ST" : "Sao Tome and Principe", "SV" : "El Salvador", "SX" : "Sint Maarten (Dutch part)", "SY" : "Syrian Arab Republic", "SZ" : "Swaziland", "TC" : "Turks and Caicos Islands", "TD" : "Chad", "TF" : "French Southern Territories", "TG" : "Togo", "TH" : "Thailand", "TJ" : "Tajikistan", "TK" : "Tokelau", "TL" : "Timor-Leste", "TM" : "Turkmenistan", "TN" : "Tunisia", "TO" : "Tonga", "TR" : "Turkey", "TT" : "Trinidad and Tobago", "TV" : "Tuvalu", "TW" : "Taiwan, Province of China", "TZ" : "Tanzania, United Republic of", "UA" : "Ukraine", "UG" : "Uganda", "UM" : "United States Minor Outlying Islands", "US" : "United States", "UY" : "Uruguay", "UZ" : "Uzbekistan", "VA" : "Holy See (Vatican City State)", "VC" : "Saint Vincent and the Grenadines", "VE" : "Venezuela, Bolivarian Republic of", "VG" : "Virgin Islands, British", "VI" : "Virgin Islands, U.S.", "VN" : "Viet Nam", "VU" : "Vanuatu", "WF" : "Wallis and Futuna", "WS" : "Samoa", "YE" : "Yemen", "YT" : "Mayotte", "ZA" : "South Africa", "ZM" : "Zambia", "ZW" : "Zimbabwe"}

#-------------------------------------------------------------------------------

def cluster_attributes (cluster):
	by_country = {}
	by_subtype = {}
	
	for seq in cluster:
		bits = seq.split ('_')
		if bits[1] in country_codes:
			hash = country_codes[bits[1]]
		else:
			hash = "Unknown"
		if hash in by_country:
			by_country [hash] += 1
		else:
			by_country [hash] = 1
			
		if bits[0] in by_subtype:
			by_subtype[bits[0]] += 1
		else:
			by_subtype[bits[0]] = 1
	
	#by_country.sort()
	#by_subtype.sort()
	
	return [by_country, by_subtype]
	
	
#-------------------------------------------------------------------------------

def report_cluster (node_count, cluster_info):
	output = cStringIO.StringIO()
	output.write ("\tCountry information:\n")
	for country in cluster_info[0]:
		output.write ("\t\t%s : %d (%g %%)\n" % (country, (cluster_info[0])[country], 100. * (cluster_info[0])[country] / node_count))
	output.write("\tSubtype information:\n")
	for subtype in cluster_info[1]:
		output.write ("\t\t%s : %d (%g %%)\n" % (subtype, (cluster_info[1])[subtype], 100. * (cluster_info[1])[subtype] / node_count))
	formatted = output.getvalue()
	output.close()
	return formatted

#-------------------------------------------------------------------------------


arguments = argparse.ArgumentParser(description='Read filenames.')

arguments.add_argument('-i', '--input',   help = 'Input CSV file with cluster definitions. Must be a CSV file with two columns: ID1,cluster ID (integer).', required = True)
arguments.add_argument('-t', '--threshold', help = 'Only report clusters with at least this many sequences')
arguments.add_argument('-a', '--alignment',   help = 'An alignment file with corresponding sequences (must be readable by HyPhy)')
arguments.add_argument('-x', '--extract',   help = 'Extract information about the cluster with a given numeric ID')
arguments.add_argument('-o', '--output',   help = 'Output sequences grouped by cluster to this directory, naming them seqid.nex (default is to output to the current directory)')

settings = arguments.parse_args()

try:
	settings.input = open (settings.input, 'r')
except IOError:
	print "Failed to open '%s' for reading" % (settings.input)
	raise

	
if settings.alignment != None:
	try:
		test = open (settings.alignment, 'r')
		test.close() 
	except IOError:
		print "Failed to open '%s' for reading" % (settings.alignment)
		raise
		
if settings.output != None:
	if os.path.isdir (settings.output) == False:
		raise Exception ("'%s' is not an existing directory " % (settings.output))
else:
	settings.output = os.getcwd()


if settings.threshold == None:
	settings.threshold = 1
else:
	settings.threshold = int (settings.threshold)
	
if settings.extract != None:
	settings.extract = int (settings.extract)

clusterReader = csv.reader(settings.input)
header = clusterReader.next ()

if len (header) != 2:
	raise Exception ('Expected a .csv file with 2 columns as input')
	
clusters = {}

for line in clusterReader:
	clusterID = int(line[1])
	if settings.extract == None or clusterID == settings.extract:
		if clusterID not in clusters:
			clusters[clusterID] = [line[0],]
		else:
			clusters[clusterID].append (line[0])
		
write_these = {}	
cluster_info = {}
		
for cluster in clusters:
	if len (clusters[cluster]) >= settings.threshold:
		print "\nCluster %d with %d sequences" % (cluster, len(clusters[cluster]))
		info = report_cluster(len(clusters[cluster]), cluster_attributes (clusters[cluster]))
		print info
		write_these [cluster] = ";".join(clusters[cluster])
		cluster_info [cluster] = info
		
if settings.alignment != None:
	print "Processing the alignment through HyPhy"
	hy_instance = hy.HyphyInterface ();
	
	script_path = os.path.realpath(__file__)
	hbl_path =  os.path.join(os.path.dirname(script_path), "HBL", "ExtractACluster.bf")
	
	hy_instance.queuevar ('alignment',  os.path.realpath(settings.alignment))
	hy_instance.queuevar ('cluster_definitions',  write_these)
	hy_instance.queuevar ('output_path', os.path.realpath(settings.output))
	hy_instance.runqueue (batchfile = hbl_path)
	print hy_instance.stdout

	
