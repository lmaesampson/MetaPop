from multiprocessing import Process
import time
from subprocess import call
import numpy as np
import itertools
import os
import zipfile

# Must be a more pythonic way:
def callmp(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25):
    call([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25])


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

    argdict = { 'beta'            : 16,
                'maxvaccprob'     : 0.8,
                'minvaccprob'     : 0.,   #this serves as the standard deviation in max. vacc. prob.
                'timesteps'       : 2000,
                'patchpop'        : 50000,
                'birthrate'       : 10.,
                'popstddev'       : 0.,
                'images'          : 0,
                'logs'            : 1,
                'fn'              : "rand_beta_smallgrid",
                'iter'            : 1,
                'stochastic'      : 0
                }

    run_process_limited_time(2000,argdict)



if __name__ == "__main__":
    main()
