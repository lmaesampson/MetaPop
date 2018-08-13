from multiprocessing import Process, Queue
import time
from subprocess import call
from subprocess import check_output
import numpy as np
import itertools
import os
import zipfile
import math



# Must be a more pythonic way:
def callmp(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,done):
    res = check_output([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29])
    done.put(res)

def run_process_limited_time(max_seconds,argdict,done):

    a = ["./metapop"]
    for key, value in argdict.items():
        a.append(key)
        a.append(str(value))

    print(a)
    a.append(done)

    # init Process
    p = Process(target=callmp, args=(a[0:(len(a))])) # 28

    # start process
    p.start()

    run_time = 0
    sleep_time = 1

    while True:
        time.sleep(sleep_time)
        run_time += sleep_time
        if not p.is_alive():
            break

        if run_time > max_seconds:
            print("Process was terminated, due to exceeding max time")
            p.terminate()
            break


def main():

    # hyperparameters
    npop = 8 # population size
    nvars =  3
    sigma = 0.01 # 0.01 noise standard deviation
    v_sigma = np.array([
        0.01,
        0.01,
        0.01
        ])
    alpha = 0.001 # 0.001learning rate

    # initialize memory for a population of w's, and their rewards
    N = np.random.randn(npop, nvars) # samples from a normal distribution N(0,1)


    argdict = { 'beta'            : 2.1,
                'maxvaccprob'     : 0.3,
                'minvaccprob'     : 0.0,
                'timesteps'           : 120,
                'patchpop'        : 10000,
                'birthrate'       : 40.0,
                'popstddev'       : 0.001,
                'images'          : 0,
                'logs'            : 0,
                'conprob'     : 0.01,  #0.08
                'constren'     : 0.01,  #0.086
                'fn'              : "payload",
                'iter'            : 5,
                'stochastic'      : 1
                }


    trns = 0

    # start the optimization
    w = np.array([
        argdict['beta'],
       # argdict['conprob'],
        argdict['maxvaccprob'],
        argdict['popstddev']
    ])

    done = Queue()
    #M = np.array([[0.001,0.001,0.001,0.001],[0.001,0.001,0.001,0.001],[1,1,1,1],[1,1,1,1]])

    # Should extend this to make better use of multiprocess
    #while trns < len(mparams):
    while trns in range(1000000):
        #print(trns)

        redflag = 0
        prvN = N
        N = np.random.randn(npop, nvars) # samples from a normal distribution N(0,1)

        R = np.zeros(npop)
        for j in range(npop):
            w_try = w + v_sigma*N[j] # jitter w using gaussian of sigma 0.1
            #print(j)
            #w_try = abs(w_try)
            for a in range(len(w_try)):
                # if w_try[a] > 1.0:
                #     w_try[a] = 1.0
                if w_try[a] < 0.0:
                    w_try[a] = 0.0

            print("w_try:",w_try)

            #argdict['constren']    = w_try[0]
            argdict['beta']   = w_try[0]
            argdict['maxvaccprob']  = w_try[1]
            argdict['popstddev'] = w_try[2]
            argdict['fn']          = "afn" + str(trns)
            argdict['logs'] = 0
            #if(j==npop-1):
            #    argdict['logs'] = 1

            #print(argdict)
           # print(list(argdict.values()))

            run_process_limited_time(1200,argdict,done)
            last_done = done.get()
            # handle reward here:
            pres = str(last_done)
            print(pres[2:len(pres)-1])
            tokens = str(pres).split(' ')

            # tksd = tokens[len(tokens)-2]
            # ftksd = float(tksd)  # if iterates, will be std dev. Otherwise, case count
            # tkm = tokens[len(tokens)-3]
            # ftkm = float(tkm) # if iterates, will be mean.
            # tkmax = tokens[len(tokens)-4]
            # ftkmax = float(tkmax) # if iterates, will be mean.
            # tkmin = tokens[len(tokens)-5]
            # ftkmin = float(tkmin) # if iterates, will be mean.
            # tsck = tokens[len(tokens)-6]
            # ftsck = float(tsck) # if iterates, will be mean.


           # reward = (ftkmax - ftkmin)/ftkmax
            reward = float(tokens[len(tokens)-2])
            #reward = (1.0 + ftksd) * (1.0 + ftsck) - 1.0
            print("reward $reward",reward)
            #reward = -ftkm

        # if(ftkmax>0.0):
        #     ltok = (ftkmax - ftkmin)/(ftkmax)
        #     print("reward:",ltok)

            if(reward == 0.0 or math.isinf(reward)):
                redflag = 1
                reward = np.random.normal(0.0, 0.001, 1)
            R[j] = reward
            #print(R[j])

          # standardize the rewards to have a gaussian distribution
        print("R:",R)
        A = np.random.randn(npop)
        print("A:",A)
        if((np.std(R)>0.0)==False):
            redflag = 1
        if((np.std(R)>0.0)==True):
            A = (R - np.mean(R)) / np.std(R)
        #print("A:",A)
        # perform the parameter update. The matrix multiply below
        # is just an efficient way to sum up all the rows of the noise matrix N,
        # where each row N[j] is weighted by A[j]
        print("A:",A)

        if(redflag == 1):
            print("redflag")

        w = w + alpha/(npop*sigma) * np.dot(N.T, A)
        # need a queue of old w to get out of trouble
        print("w:",w)

        trns = trns + 1
                    # add system calls to ffmpeg

if __name__ == "__main__":
    main()


#     {'beta': 7.9938039771812974, 'maxvaccprob': 0.36938294758894308, 'minvaccprob': 0.0, 'timesteps': 10, 'patchpop': 11744.457220245346, 'birthrate': 0.6834794522665
# 4827, 'popstddev': 3.6, 'images': 0, 'logs': 0, 'fn': 'afn9', 'iter': 100, 'stochastic': 1}
# ['./metapop', 'beta', '7.99380397718', 'maxvaccprob', '0.369382947589', 'minvaccprob', '0.0', 'timesteps', '10', 'patchpop', '11744.4572202', 'birthrate', '0.6834
# 79452267', 'popstddev', '3.6', 'images', '0', 'logs', '0', 'fn', 'afn9', 'iter', '100', 'stochastic', '1']


#./metapop beta 7.99380397718 maxvaccprob 0.369382947589 minvaccprob 0.0 timesteps 10 patchpop 11744.4572202 birthrate 0.6834 popstddev 3.6 images 0 logs 0 fn afn9 iter 100 stochastic 1
