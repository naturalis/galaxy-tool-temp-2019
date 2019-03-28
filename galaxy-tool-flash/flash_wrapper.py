#!/usr/bin/python
"""
FLASH   V1.0    martenhoogeveen@naturalis.nl

This is a wrapper for the tool FLASH (https://ccb.jhu.edu/software/FLASH/).
FLASH can merge paired-end sequencing data. The wrapper handles a zip file with one or more pairs.
"""
import sys, os, argparse
import glob
import string
from Bio import SeqIO
from subprocess import call, Popen, PIPE

# Retrieve the commandline arguments
parser = argparse.ArgumentParser(description='')
requiredArguments = parser.add_argument_group('required arguments')
#Inputfile in zip format
requiredArguments.add_argument('-i', '--input', metavar='input zipfile', dest='inzip', type=str,
                               help='Inputfile in zip format', required=True)
#output folder
requiredArguments.add_argument('-of', '--folder_output', metavar='folder output', dest='out_folder', type=str,
                               help='Folder name for the output files', required=True)
#A zipfile can contain FASTQ or gzip files,the user needs to set this parameter.
requiredArguments.add_argument('-t', '--input_type', metavar='FASTQ or GZ input', dest='input_type', type=str,
                               help='Sets the input type, gz or FASTQ', required=True)
#This parameter determines what happens with the reads that do not merge. This are the three options:
#discard: non-merged reads are discarded
#add: forward reads are added to the fastq file with merged reads
#seperate: the merged and forward non-merged reads will be outputted in seprated files
requiredArguments.add_argument('-f', '--forward', metavar='Add non-megerd forward reads', dest='forward', type=str,
                               help='Adds the forward reads to the merged reads file. Option can be discard, add or seperate', required=True)
#minimum bases that need to overlap to merge the reads
requiredArguments.add_argument('-m', '--min-overlap', metavar='minimum overlap', dest='minforward', type=str,
                               help='minimum overlap', required=True)
#The allowed mismatch percentage for merging
requiredArguments.add_argument('-x', '--mis-ratio', metavar='mismatch ratio', dest='mismatch', type=str,
                               help='mismatch ratio', required=True)
#maximum number of bases that are allowed to overlap
requiredArguments.add_argument('-M', '--max-overlap', metavar='max overlap', dest='maxoverlap', type=str,
                               help='maxoverlap', required=True)
args = parser.parse_args()

def admin_log(tempdir, out=None, error=None, function=""):
    """
    A log file will be made and log data will be written to that file. Most of the time this is the stdout and stderror
    of the shell. In the log it says if the message in is coming from stdout or stderror.
    :param tempdir: the tempdir path that contains the log file
    :param out: stdout or out message
    :param error: stderror or error message
    :param function: name of the function or step that generated the message
    """
    with open(tempdir + "/log.log", 'a') as adminlogfile:
        seperation = 60 * "="
        if out:
            adminlogfile.write("out "+ function + " \n" + seperation + "\n" + out + "\n\n")
        if error:
            adminlogfile.write("error " + function + "\n" + seperation + "\n" + error + "\n\n")

def make_output_folders(tempdir):
    """
    Output en work folders are created. The wrapper uses these folders to save the files that are used between steps.
    :param tempdir: tempdir path
    """
    call(["mkdir","-p", tempdir])
    call(["mkdir", tempdir + "/paired_files"])
    call(["mkdir", tempdir + "/merged_files"])
    call(["mkdir", tempdir + "/output"])

def get_files(tempdir):
    """
    This function finds the file pairs. First it looks at the beginning of the filename before the seperator (R1)
    after that it looks for a file with that same beginning. The are stored in a list, and all the lists are stored
    in a dictionairy. Example: {seqfile:[seqfile1_R1_miseq.fastq,seqfile1_R2_miseq.fastq]}
    :param tempdir: the tempdir path that contains the log file
    :return: A dictionairy where the keys are a part of the file name en the values a list with the forward and reverse
    filename.
    """
    filetype = tempdir+"/paired_files/*.fastq"
    gzfiles = [os.path.basename(x) for x in sorted(glob.glob(filetype))]
    reverse=[]
    pairs={}
    for x in gzfiles:
        if x not in reverse:
            sample = x.partition("R1")[0]
            pairlist=[]
            for y in gzfiles:
                if sample == y[:len(sample)]:
                    pairlist.append(y)
                    if y[:(len(sample)+2)] == sample+"R2":
                        reverse.append(y)
            pairs[sample[:-1]] = pairlist
    return pairs

def gunzip(tempdir):
    """
    If the input zip file contains gzip files they need to be gunzipped. The files are gunzipped and placed in the
    paired files folder. The characters dash dot and space are replaced by an underscore.
    :param tempdir: tempdir path
    """
    filetype = tempdir + "/paired_files/*.gz"
    gzfiles = [os.path.basename(x) for x in sorted(glob.glob(filetype))]
    for x in gzfiles:
        call(["gunzip", tempdir + "/paired_files/" + x])
        gunzip_filename = os.path.splitext(x[:-3])
        if x[:-3] != gunzip_filename[0].translate((string.maketrans("-. ", "___"))) + gunzip_filename[1]:
            call(["mv", tempdir + "/paired_files/" + x[:-3], tempdir + "/paired_files/" +gunzip_filename[0].translate((string.maketrans("-. " , "___")))+gunzip_filename[1]])

def flash(pairs, tempdir):
    """
    This method execute the FLASH merger. The output files are written to the folder merged_files. All the output files
    will be moved to the output folder depending on the forward parameter. Later the output folder will be zipped
    :param pairs: dictionairy with the filenames.
    :param tempdir: tempdir path
    """
    for x in pairs:
        basename = pairs[x][0].split("_R1")[0]
        out, error = Popen(["flash", tempdir+"/paired_files/"+pairs[x][0], tempdir+"/paired_files/"+pairs[x][1],"-x", args.mismatch ,"-m", args.minforward, "-M", args.maxoverlap, "-d", tempdir+"/merged_files/" ,"-o", basename], stdout=PIPE, stderr=PIPE).communicate()
        admin_log(tempdir, out=out, error=error, function="flash")
        if args.forward == "add":
            with open(tempdir+"/output/"+basename+"_merged_forward.fastq", 'a') as outfile:
                call(["cat", tempdir+"/merged_files/"+basename+".extendedFrags.fastq", tempdir+"/merged_files/"+basename+".notCombined_1.fastq"],stdout=outfile)
            call(["rm", tempdir + "/merged_files/" + basename + ".notCombined_1.fastq", tempdir + "/merged_files/" + basename + ".extendedFrags.fastq"])
        if args.forward == "seperate":
            call(["mv",  tempdir+"/merged_files/"+basename+".extendedFrags.fastq", tempdir+"/output/"+basename+"_merged.fastq"])
            call(["mv",  tempdir+"/merged_files/"+basename + ".notCombined_1.fastq", tempdir + "/output/" + basename + "_forward.fastq"])
        if args.forward == "discard":
            call(["mv",  tempdir+"/merged_files/"+basename + ".extendedFrags.fastq", tempdir + "/output/" + basename + "_merged.fastq"])
            call(["rm", tempdir+"/merged_files/"+basename + ".notCombined_1.fastq"])
        call(["rm", tempdir+"/paired_files/"+pairs[x][0], tempdir+"/paired_files/"+pairs[x][1], tempdir + "/merged_files/" + basename + ".histogram", tempdir+"/merged_files/"+basename+".hist", tempdir+"/merged_files/"+basename+".notCombined_2.fastq"])

def zip_it_up(tempdir):
    """
    Files in the outfolder will be zipped and moved to the galaxy output path.
    :param tempdir: tempdir path
    :return:
    """
    call(["zip","-r","-j", tempdir+"/merged.zip", tempdir+"/output/"],stdout=open(os.devnull, 'wb'))
    call(["rm", "-rf", tempdir + "/output", tempdir + "/paired_files", tempdir + "/merged_files"])

def main():
    tempdir = args.out_folder
    make_output_folders(tempdir)
    zip_out, zip_error = Popen(["unzip", args.inzip, "-d", tempdir.strip() + "/paired_files"], stdout=PIPE,stderr=PIPE).communicate()
    admin_log(tempdir, zip_out, zip_error)
    if args.input_type == "gz":
        gunzip(tempdir)
    pairs = get_files(tempdir)
    flash(pairs, tempdir)
    zip_it_up(tempdir)

if __name__ == '__main__':
    main()



