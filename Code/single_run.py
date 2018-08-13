from multiprocessing import Process
import time
from subprocess import call
import numpy as np
import itertools
import os
import zipfile

# Must be a more pythonic way:
def callmp(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29):
    call([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29])


def run_process_limited_time(max_seconds,argdict):

    a = ["./metapop"]
    for key, value in argdict.items():
        a.append(key)
        a.append(str(value))

    print(a)

    # init Process
    p = Process(target=callmp, args=(a[0:len(a)]))

    # start process
    p.start()

    run_time = 0
    sleep_time = 1

    while True:
        time.sleep(sleep_time)
        run_time += sleep_time
        #print(run_time)
        if not p.is_alive():
            break

        if run_time > max_seconds:
            print("Process was terminated, due to exceeding max time")
            p.terminate()
            break


def main():

    argdict = { 'beta'        : 2.0,
                'maxvaccprob' : 0.5,
                'minvaccprob' : 0.1,
                'timesteps'   : 100,
                'patchpop'    : 10000,
                'birthrate'   : 40.0,
                'popstddev'   : 0.1,
                'images'      : 1,
                'logs'        : 1,
                'conprob'     : 0.001,  #0.086
                'constren'     : 0.01,  #0.086
                'fn'          : "perm",
                'iter'        : 1,
                'stochastic'  : 1
                }

    trns = 0

    # Should extend this to make better use of multiprocess
    #while trns < len(mparams):
    #while trns in range(len(mparams)):
    while trns in range(1):
        #print(trns)
        run_process_limited_time(600000,argdict)
        if(argdict['images'] == 1):
            monfs = "omontage" + str(trns) + ".png"
            mtchs = "./out_" + argdict['fn'] + "*.png"
            call(["montage", "-mode", "concatenate", mtchs, monfs])
            call(["rm", mtchs])
        if(argdict['logs'] == 1):
            filestem = "./cases_" + argdict['fn'] + "_" + str(0) # not handling iterates yet...
            filetc = filestem + ".csv" # cases_out_afn0_5.csv # must handle iterates
            filezip = filestem + ".zip" # cases_out_afn0_5.csv # must handle iterates
           # zip_name = zipfile.ZipFile(filezip, 'a')
            #zip_name.write(filetc, os.path.basename(filetc), zipfile.ZIP_DEFLATED)
          #  call(["rm", filetc])
        trns = trns + 1
        # add system calls to ffmpeg

if __name__ == "__main__":
    main()
