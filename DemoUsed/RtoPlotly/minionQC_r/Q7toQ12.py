import argparse
import multiprocessing as mp
import os
import subprocess as sup
import sys
#-------------
"""
    Author: Gray
    Date:   2018-10-02
    Aim:    Run MinIONQC.r at -q 7 to 12
"""


def JOB(workTuple):
    Q, inputSummarySeq, workDir, Rscript, FileMarker = workTuple
    os.makedirs(os.path.join(workDir,"Q_"+str(Q)+"out"))

    if os.name == "nt":  #Windows
        EXE = "Rscript.exe"
    elif os.name == "posix": #Unix
        EXE = "Rscript"
    else:
        print("unknow OS")
        os._exit()

    cmd1 = [EXE, Rscript, "-i", inputSummarySeq,
    "-o", os.path.join(workDir,"Q_"+str(Q)+"out"), "-q", str(Q),
    "-m", FileMarker]
    sup.call(cmd1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=sys.argv[0]+
                                     ": Run MinIONQC.R in Q7 to Q12")
    parser.add_argument('--script', dest="MainScript",
                    help="The MinIONQC script, for different modified version",
                    default=r"D:\ProjectWork\20181002_ONT_MinION_QC_modifying\MinIONQC_modByGray_20181002.r")
    parser.add_argument('-i', '--input', dest="INPUT", required=True,
                    help="input inputSummaryText from albacore")
    parser.add_argument('-m', '--marker', dest="Marker", required=True,
                    help="The marker as tag for yaml result sheets")
    parser.add_argument('-o', '--outdir', dest="OUTDIR", required=True,
                    help="output work dir, auto-create Q7 to Q12 folders")
    parser.add_argument("--core", dest="CORE", default=1)


    results = parser.parse_args()
    INPUT = results.INPUT
    WORKDIR = results.OUTDIR
    Rscript = results.MainScript
    FileMarker = results.Marker
    core = int(results.CORE)
    mp.freeze_support()
    worklist = [ (Q, INPUT, WORKDIR, Rscript, FileMarker) for Q in range(7,13)]
    pool = mp.Pool(core) #os.cpu_count()
    pool.map(JOB, worklist)
