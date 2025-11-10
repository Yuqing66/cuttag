import linecache
import argparse
import os.path

# Argument parsing:
def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-odir', type=str, required=True,
						help='output folder path')
	parser.add_argument('-f', type=str, required=True,
						help='filename')
	parser.add_argument('-endbp', type=int, required=True,
						help='end range that count as duplicates')
	args = parser.parse_args()
	return args

# Add new dup to the dictionary for storing duplicates info
def adddupdict(dupdict, umi):
	info = umi.split(':')
	dupdict['start'].add(info[1])
	dupdict['end'].add(int(info[2]))
	dupdict['ab1'].add(info[3])
	dupdict['ab2'].add(info[4])

# is this line duplicates?
def isdup(dupdict,umi,endbp):
	decision = False
	info = umi.split(":")
	start = info[1]
	end = int(info[2])
	ab1 = info[3]
	ab2 = info[4]
	if ab1 in dupdict['ab1'] and ab2 in dupdict['ab2'] and start in dupdict['start']:
		if end in dupdict['end']:
			decision = True
		else:
			endextent = []
			for i in range(1,endbp+1):
				endextent += [x+i for x in dupdict['end']]
				endextent += [x-i for x in dupdict['end']]
			if end in endextent:
				decision = True
	return decision

def outFormat(dups):
	locab = {}
	cb = {}
	for umi, count in dups.items():
		info = umi.split(':')
		umiL = '%s:%s:%s:%s:%s' % (info[0],info[1],info[2],info[3],info[4])
		umiR = '%s' % (info[5])
		if locab.has_key(umiL):
			locab[umiL] += count
		else:
			locab[umiL] = count
		if cb.has_key(umiR):
			cb[umiR] += count
		else:
			cb[umiR] = count
	locab_most = max(locab, key=locab.get)
	record = '%s\t%s\t%s\n' % (locab_most, ':'.join(cb.keys()), ':'.join(map(str,cb.values())))
	return record

def run(args):
	odir = args.odir
	inFile = args.f
	inf = open(inFile, 'r')
	endbp = args.endbp

	# antibody barcodes frequency within each dups
	freqFile = os.path.join(odir, inFile.replace('.txt','_freq.txt'))
	freq = open(freqFile, 'wb')

	dup_bc_counter = 0

	# dedup
	dupdict = {'start':set(),'end':set(),'ab1':set(),'ab2':set()}
	dups = {}
	# first line
	line = inf.readline().split()
	count = int(line[0])
	umi = line[1]
	adddupdict(dupdict,umi)
	dups[umi] = count

	while 1:
		line = inf.readline().split()
		if not line:
			break
		count = int(line[0])
		umi = line[1]
		if isdup(dupdict, umi, endbp):
			adddupdict(dupdict,umi)
			dups[umi] = count
		else:
			# write the previous dup to file.
			record = outFormat(dups)
			freq.write(record)
			if len(record.split('\t')[2].split(':')) > 1:
				dup_bc_counter += 1

			# start a new dups
			dupdict = {'start':set(),'end':set(),'ab1':set(),'ab2':set()}
			dups = {}
			adddupdict(dupdict,umi)
			dups[umi] = count

	# write the last dups
	record = outFormat(dups)
	freq.write(record)
	if len(record.split('\t')[2].split(':')) > 1:
				dup_bc_counter += 1

	freq.close()
	
	counter = open(os.path.join(odir, 'counter.txt'), 'wb')
	counter.write('dup_bc_counter = %s' % dup_bc_counter)
	counter.close()
	return

if __name__ == "__main__":
	args = get_args()
	run(args)

