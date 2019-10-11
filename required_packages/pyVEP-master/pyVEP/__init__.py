
import urllib
import logging
import requests

__version__ = '0.0.1'

'''
Example curl command:
curl 'http://rest.ensembl.org/vep/Homo_sapiens/hgvs/9%3Ag.22125504G%3EC' \
-H 'Host: rest.ensembl.org' \
-H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.7; rv:42.0) Gecko/20100101 Firefox/42.0' \
-H 'Accept: application/json, text/javascript, */*; q=0.01' \
-H 'Accept-Language: en-US,en;q=0.5' \
--compressed \
-H 'Referer: http://asia.ensembl.org/Tools/VEP?redirect=no' \
-H 'Origin: http://asia.ensembl.org' \
-H 'Connection: keep-alive'
'''

'''
TODO:
* Randomize headers 
    Similar to this: https://github.com/galkan/tools/blob/master/others/programming/python/random-http-headers-urllib.py 
'''

class PyVEPException(Exception):
	pass


def get_variant_type(variant):
	'''
	Simple heuristics to get the variant type

	1 182712 182712 A/C 1  --> 1:182712-182712:1/C
	3 319781 319781 A/- 1  --> 3:319781-319781:1/-
	19 110748 110747 -/T 1 --> 19:110748-110747:1/T
	1 160283 471362 DUP 1 --> 1:160283-471362:1/DUP
	1 1385015 1387562 DEL 1 --> 1:1385015-1387562:1/DEL  

	1 182712 . A C . . . --> 1:182712-182712:1/C
	3 319780 . GA G . . . --> 3:319781-319781:1/- 
	3 319780 . GAA G . . . --> 3:319781-319782:1/- 
	19 110747 . G GT . . . --> 9:110748-110747:1/T 
	19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
	1 160283 sv1 . <DUP> . . SVTYPE=DUP;END=471362 . --> UNABLE TO GENERATE PREVIEW!
	1 1385015 sv2 . <DEL> . . SVTYPE=DEL;END=1387562 . --> UNABLE TO GENERATE PREVIEW!
	'''

	new_variant = variant
	variant_splitted = variant.split()
	if len(variant_splitted) == 5:
		# Ensembl default

		ref_alt = variant_splitted[3].split('/')
		if len(ref_alt) == 2:
			if ref_alt[1] == '-':
				new_variant = ''.join([
					variant_splitted[0], 
					':', 
					variant_splitted[1], 
					'-', 
					variant_splitted[2],
					':',
					variant_splitted[4],
					'/',
					'-',
					])
			else:
				new_variant = ''.join([
					variant_splitted[0], 
					':', 
					variant_splitted[1], 
					'-', 
					variant_splitted[2],
					':',
					variant_splitted[4],
					'/',
					ref_alt[1],
					])
		elif variant_splitted[3] in ['DUP', 'DEL']:
			new_variant = ''.join([
				variant_splitted[0], 
				':', 
				variant_splitted[1], 
				'-', 
				variant_splitted[2],
				':',
				variant_splitted[4],
				'/',
				variant_splitted[3],
				])
		else:
			message = '%s notation is not supported..' % (variant_splitted[3])
			logging.error(message)
			raise PyVEPException(message)

		return 'region', new_variant

	if len(variant_splitted) == 8:
		if len(variant_splitted[3]) == 1 and len(variant_splitted[4]) == 1:
			#Simple SNP
			new_variant = ''.join([
				variant_splitted[0],
				':',
				variant_splitted[1],
				'-',
				variant_splitted[1],
				':',
				'1', # Not sure what this is..
				'/',
				variant_splitted[4],
				])
		elif len(variant_splitted[3]) > len(variant_splitted[4]):
			# Deletion
			new_variant = ''.join([
				variant_splitted[0],
				':',
				str(int(variant_splitted[1]) + len(variant_splitted[4])),
				'-',
				str(int(variant_splitted[1]) + len(variant_splitted[3]) - len(variant_splitted[4])),
				':',
				'1', # Not sure what this is..
				'/',
				'-',
				])
		elif len(variant_splitted[3]) < len(variant_splitted[4]):
			# Insertion
			# 19 110747 . G GT . . . --> 9:110748-110747:1/T 
			# 19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
			new_variant = ''.join([
				variant_splitted[0],
				':',
				str(int(variant_splitted[1]) + len(variant_splitted[3])),
				'-',
				variant_splitted[1],
				':',
				'1', # Not sure what this is..
				'/',
				variant_splitted[4][len(variant_splitted[3]):],
				])
		else:
			raise PyVEPException('Could not parse VCF variant: %s' % (variant))

		return 'region', new_variant

	if variant[0:2] == 'rs':
		# Variant identifier
		return 'id', new_variant

	if ':' in variant:
		# HGVS notation 
		return 'hgvs', new_variant

	if len(variant_splitted) > 1:
		logging.error('Could not recognize variant: %s' % (variant))

	return 'id', new_variant # Hoping for the best

def VEP(variant, assembly, variant_type=None):
	
	url_pattern_grch38 = 'http://rest.ensembl.org/vep/Homo_sapiens/{variant_type}/{variant}'
	url_pattern_grch37 = 'http://grch37.rest.ensembl.org/vep/Homo_sapiens/{variant_type}/{variant}'

	preprocess_assembly = lambda x : x.split('.')[0].lower()

	if preprocess_assembly(assembly) in ['grch38', 'hg38']:
		url_pattern = url_pattern_grch38
		header_version = 'rest'
	elif preprocess_assembly(assembly) in ['grch37', 'hg19']:
		url_pattern = url_pattern_grch37
		header_version = 'grch37.rest'
	else:
		message = 'Unknown assembly: %s Accepted values: grch38, hg38, grch37, hg19' % (str(assembly))
		logging.error(message)
		raise ValueError(message)

	headers = {
		'Host': '%s.ensembl.org' % (header_version),
		'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.7; rv:42.0) Gecko/20100101 Firefox/42.0',
		'Accept': 'application/json, text/javascript, */*; q=0.01',
		'Accept-Language': 'en-US,en;q=0.5',
		'Referer': 'http://asia.ensembl.org/Tools/VEP?redirect=no',
		'Origin': 'http://asia.ensembl.org',
		'Connection': 'keep-alive',
	}


	if variant_type is None:
		variant_type, variant = get_variant_type(variant)
	if not type(variant_type) is str:
		raise PyVEPException('variant_type should be str')
	if not type(variant).__name__ in ['str', 'unicode']:
		raise PyVEPException('variant should be str or unicode')

	if str(variant_type) not in ['region', 'id', 'hgvs']:
		raise ValueError('Unknown variant_type value: %s . Accepted values: region, id, hgvs' % (str(variant_type)))

	logging.debug('VAP Variant: %s' % (variant))
	variant_quoted = urllib.quote(variant)
	url = url_pattern.format(
		variant=variant_quoted, variant_type=variant_type)

	logging.debug('URL: %s' % (url))
	r = requests.get(url, headers=headers)
	json_text = r.json()

	return json_text

def test():
	import json
	
	logging.basicConfig(level=logging.DEBUG)

	def try_variant(variant):
		logging.info('Trying variant: %s' % (variant))
		r = VEP(variant, 'grch38')
		logging.info('Result: %s' % json.dumps(r, indent=4))

	testing_variants = {
		'ENSEMBL default variant': [
			'1 182712 182712 A/C 1',
			'2 265023 265023 C/T 1',
			'3 319781 319781 A/- 1',
			'19 110748 110747 -/T 1',
			'1 160283 471362 DUP 1',
			'1 1385015 1387562 DEL 1',
		],
		'VCF': [
			'1 182712 . A C . . .',
			'3 319780 . GA G . . .',
			'3 319780 . GAA G . . .',
			'19 110747 . G GT . . .',
			'19 110747 . G GTT . . .',
#			'1 160283 sv1 . <DUP> . . SVTYPE=DUP;END=471362 .',
#			'1 1385015 sv2 . <DEL> . . SVTYPE=DEL;END=1387562 .',
		],
		'Variant identifiers': [
			'rs699',
			'rs144678492',
			'COSM354157',
		],
		'HGVS notations': [
			'AGT:c.803T>C',
			'9:g.22125504G>C',
			'ENST00000003084:c.1431_1433delTTC',
			'19:g.110747_110748insT',
			'LRG_101t1:c.1019T>C',
		],
	}

	for variant_type in testing_variants:
		logging.info('TESTING TYPE: %s' % (variant_type))
		for variant in testing_variants[variant_type]:
			try_variant(variant)
