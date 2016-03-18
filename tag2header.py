#!/usr/bin/python


from optparse import OptionParser
import sys
import os


import sys
import os
import time
import logging
import gzip
import random
import Bio.SeqIO
import dms_tools
import dms_tools.parsearguments
import dms_tools.file_io
import dms_tools.utils
import dms_parse_modified
import dms_utils_modified
import Levenshtein as L

usage = "usage: %prog [options] \n"

prog="this prog"

parser = OptionParser(usage=usage)

parser.add_option("--outprefix",help="outdir/outputfile prefix")
parser.add_option("--r1",help="R1 file fastq.gz")
parser.add_option("--r2",help="R2 file fastq.gz")


parser.add_option("--minq",default=20,help="reads min quality")
parser.add_option("--minqbarcode",default=20,help="barcode min quality")
parser.add_option("--maxlowqfrac",default="0.025",help="max low quality base fraction")
parser.add_option("--barcodelength",default="15",help="barcode length")
parser.add_option("--mergebarcode",default="0",help="whether merge barcod[0]")
parser.add_option("--barcodeerror",default="1",help="barcode error number")

(opts,args)=parser.parse_args()


def hdr_rename_fxn(read_title, read_tag):
    # put barcode to the read name
    illumina = read_title.split(" ")[0].split(":")

    if len(illumina) == 7:
        #Illumina CASAVA >=1.8
        #e.g. @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACGcd
        return "%s|%s" % (read_title.split(" ")[0], read_tag)
    elif len(illumina) == 5:
        #Illumina CASAVA >=1.4?
        #e.g. @HWUSI-EAS100R:6:73:941:1973#ATCGAT/1
        read_title = read_title.replace(' ', '_')
        return "%s|%s%s/%s" % (read_title.split('/')[0], read_tag, read_title.split('/')[1])
    else :
        raise ValueError("Unknown read name format: %s" % read_title)


def main():
    # Parse command line arguments
    (opts,args)=parser.parse_args()
    r1file=[opts.r1]
    r2file=[opts.r2]
    minq=int(opts.minq)
    minqbarcode=int(opts.minqbarcode)
    maxlowqfrac=float(opts.maxlowqfrac)
    barcodelength=int(opts.barcodelength)
    mergebarcode=int(opts.mergebarcode)
    barcodeerror=int(opts.barcodeerror)
    outtag=opts.outprefix
    # set up logging

    logging.shutdown()
    logfile = "%s.log" % outtag

    if os.path.isfile(logfile):
        os.remove(logfile)
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
    logger = logging.getLogger(prog)
    logfile_handler = logging.FileHandler(logfile)
    logger.addHandler(logfile_handler)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logfile_handler.setFormatter(formatter)
    logger.info("Beginning execution of %s in directory %s\n" % (prog, os.getcwd()))
    logger.info("Progress will be logged to %s" % logfile)

    # define file names and delete existing files
    r1name = '%s.R1.fq.gz' %outtag
    r2name = '%s.R2.fq.gz' %outtag
    r2wibar = '%s.Rwithbar.fq.gz' %outtag
    outR1 = gzip.open(r1name, 'w')
    outR2 = gzip.open(r2name, 'w')
#    outR2wbar=gzip.open(r2wibar, 'w')
    
    outfiles = [outR1,outR2]
    try:		
        # check on read files

        # collect reads by barcode while iterating over reads
        readcategories = ['total read pairs', 'read pairs that fail Illumina filter', 'low quality read pairs']
        n = dict([(category, 0) for category in readcategories])
        barcodes = {}
        gzipped=True
        print r1file
        for read_tup in dms_tools.file_io.IteratePairedFASTQ(r1file, r2file, gzipped, applyfilter=False):
            n['total read pairs'] += 1
            if read_tup:
                (name, r1, r2, q1, q2) = read_tup
        #        print name
                qcheckedreads = dms_utils_modified.CheckReadQuality(r1, r2, q1, q2, minq, minqbarcode, maxlowqfrac, barcodelength)
                if qcheckedreads:
                    (r1, q1, r2, q2) = qcheckedreads
                    barcode =  r2[ : barcodelength]

                    newname = hdr_rename_fxn(name, barcode)
                    outR1.write('@%s\n%s\n+\n%s\n' % (newname,r1,q1))
                    outR2.write('@%s\n%s\n+\n%s\n' % (newname,r2,q2))
                    
                    if barcode in barcodes:
                        barcodes[barcode].append((name, r1, q1, r2, q2))
                    else:
                        if mergebarcode:
                            flag=0
                            for existbarcode in barcodes.iterkeys():
                                if L.distance(existbarcode,barcode) <= barcodeerror:
                                    barcodes[existbarcode].append((name, r1, q1, r2, q2))
                                    flag=1
                                    break
                            if flag==0:
                                barcodes[barcode] = [(name, r1, q1, r2, q2)]

                        else:
                            barcodes[barcode] = [(name, r1, q1, r2, q2)]
                else:
                    n['low quality read pairs'] += 1
            else:
                n['read pairs that fail Illumina filter'] += 1
            if n['total read pairs'] % 1e5 == 0:
                logger.info('Reads parsed so far: %d' % n['total read pairs'])

        logger.info('Finished parsing all %d reads; ended up with %d unique barcodes.\n' % (n['total read pairs'], len(barcodes)))

        logger.info('Now examining the %d barcodes to see if reads meet criteria for retention' % len(barcodes))


    except:
        logger.exception('Terminating %s at %s with ERROR' % (prog, time.asctime()))
        for f in outfiles:
            try:
                f.close()
            except:
                pass
        for f in outfilenames:
            if os.path.isfile(f):
                os.remove(f)
                logger.exception('Deleting file %s' % f)
    else:
        logger.info('Successful completion of %s at %s' % (prog, time.asctime()))
    finally:
        for f in outfiles:
            try:
                f.close()
            except:
                pass
        logging.shutdown()


if __name__ == '__main__':
    main() # run the script
