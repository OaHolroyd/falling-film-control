# Things that this code should be able to do:
#   - set up JSON file for given RE, CA, etc.
#   - run any of the three options and compute analytic and estimated
#     damping rates
#   - plot hstats for any combination of the three options
#   - plot damping rates accross multiple runs

import os
import copy
import math
import json
import datetime
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI


def mp_worker(inputs):
    """worker function for pool"""
    cmd = inputs[0]
    cwd = inputs[1]
    idx = inputs[2]
    N = inputs[3]
    pid = inputs[4]
    ip = inputs[5]
    print(f" [{pid}] {idx}, {ip}/{N} (started {datetime.datetime.now()})", flush=True)
    try:
        try:
            # timeout after 6 hours
            subprocess.run([cmd],
                           stderr=subprocess.DEVNULL,
                           stdout=subprocess.DEVNULL,
                           cwd=cwd, timeout=6*3600)
        except:
            pass

        # remove so that it won't rerun on restart
        subprocess.run(["rm", cmd],
                       stderr=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL,
                       cwd=cwd)

        # remove unrequired data
        subprocess.run(["rm", "-r", "dump"],
                       stderr=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL,
                       cwd=cwd)
    except FileNotFoundError:
        pass


class FilmModel():
    """
    Class to contain methods for setting up and running the C code
    """
    def __init__(self, Re=1, Ca=0.1):
        super().__init__()

        # create parameters and set to correspond to Re/Ca
        self.params = {
                           "DOMAIN": {
                             "h0": 0.000018989145744046399526063052252081,
                             "Lx": 30.0,
                             "Ly": 8.0,
                             "theta": 1.047197551,
                             "tmax": 25,
                             "t0": 0
                           },

                           "PHYSICAL": {
                             "rho_l": 10058.938842622002416886564380943,
                             "rho_g": 10.0,
                             "mu_l": 1.0e-3,
                             "mu_g": 1.0e-5,
                             "gamma": 0.00015705930063439097934322589390558,
                             "grav": 10.0
                           },

                           "SOLVER": {
                             "level": 8,
                             "dtout": 0.1,
                             "output": 0
                           },

                           "CONTROL": {
                             "M": "5",
                             "P": "5",
                             "start": 5.0,
                             "width": 0.1,
                             "alpha": 1.0,
                             "del": 0.0,
                             "mu": 0.5,
                             "rom": "benney",
                             "strategy": "lqr",
                             "exact_flux": 1
                           }
                       }
        self.set_params(Re, Ca)

        # create space for results
        self.data0 = {"analytic_be": None,
                      "analytic_wr": None,
                      "ns": None,
                      "wr": None,
                      "benney": None}
        self.data1 = {"ns": None,
                      "wr": None,
                      "benney": None}
        self.rates = {"analytic_be": None,
                      "analytic_wr": None,
                      "ns": None,
                      "wr": None,
                      "benney": None}

    def set_params(self, Re, Ca):
        """
        Set the various parameters to get a given Re/Ca pair
        """
        T = 1.20904306e-3  # time scale
        grav = self.params["PHYSICAL"]["grav"]
        sinth = math.sin(self.params["DOMAIN"]["theta"])
        mu = self.params["PHYSICAL"]["mu_l"]
        self.Re = Re
        self.Ca = Ca
        self.params["DOMAIN"]["h0"] = (Re*T*T*grav*sinth)/2.0
        self.params["PHYSICAL"]["rho_l"] = 4*mu/(Re*(T**3)*((grav*sinth)**2))
        self.params["PHYSICAL"]["rho_g"] = self.params["PHYSICAL"]["rho_l"]/1000
        self.params["PHYSICAL"]["gamma"] = (Re*T*grav*mu*sinth)/(2*Ca)

    def load_params(self, file):
        """
        Load params from a file
        """
        with open(file, 'r', encoding='utf-8') as f:
            self.params = json.load(f)

    def dump_params(self, file="params.json"):
        """
        dump the parameters to a json file ready to be run
        """
        self.params["CONTROL"]["M"] = int(self.params["CONTROL"]["M"])
        self.params["CONTROL"]["P"] = int(self.params["CONTROL"]["P"])
        with open(file, 'w', encoding='utf-8') as f:
            json.dump(self.params, f)

    def run_model(self, model=None):
        """
        Runs the code for a given model to generate the data
        """
        # write the JSON to the input file
        self.dump_params()

        cmd = ""
        if model == "ns":
            cmd = "./film-ns"
        elif model == "wr":
            cmd = "./film-wr"
        elif model == "benney":
            cmd = "./film-benney"
        else:
            cmd = "./film-benney & ./film-wr ./film-ns; wait"

        # run the correct code
        count = 0
        while count < 10:
            count += 1
            try:
                subprocess.run([cmd],
                               stderr=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                print("  failed, tring again")

    def multi_run(self, REs, CAs, rom="wr", uid=1):
        """
        Runs the code over all of the supplied parameters. The lengths must
        either match or one must have a single value
        """
        # get length and convert to lists if required
        try:
            N = len(REs)
        except TypeError:
            N = 1
            REs = [REs]
        try:
            M = len(CAs)
        except TypeError:
            M = 1
            CAs = [CAs]

        # ensure equal length
        if N != M:
            if N == 1:
                REs = 0*CAs + REs[0]
                N = M
            elif M == 1:
                CAs = 0*REs + CAs[0]
            else:
                return

        # get comms details
        comm = MPI.COMM_WORLD
        my_rank = comm.Get_rank()
        p = comm.Get_size()
        M = 0  # number of runs this task has to perform
        for i in range(N):
            if i % p != my_rank:
                continue
            M = M + 1

        # set up "queue"
        work_dir = f"py-work-{rom}-{uid:04d}"
        output_dir = f"py-out-{rom}-{uid:04d}"
        queue_dir = f"{work_dir}/queue"
        if my_rank == 0:
            os.system(f"mkdir -p {output_dir}")
            os.system(f"mkdir -p {work_dir}")
            if not os.path.isdir(f"mkdir -p {queue_dir}"):
                os.system(f"mkdir -p {queue_dir}")
                for i in range(N):
                    os.system(f"touch {queue_dir}/{i:04d}.txt")
        comm.Barrier()

        while True:
            # find next task
            ongoing = 0
            for i in range(N):
                try:
                    subprocess.run(["rm", f"{work_dir}/queue/{i:04d}.txt"],
                                   stderr=subprocess.DEVNULL,
                                   stdout=subprocess.DEVNULL,
                                   check=True)
                    ongoing = 1

                    # set up directories, executable and parameter file
                    os.system(f"mkdir {work_dir}/run-{i:04d}")
                    os.system(f"mkdir {work_dir}/run-{i:04d}/out")
                    os.system(f"mkdir {work_dir}/run-{i:04d}/dump")
                    os.system(f"cp film-ns {work_dir}/run-{i:04d}/film-ns")
                    self.params["CONTROL"]["rom"] = rom
                    self.set_params(REs[i], CAs[i])
                    self.dump_params(file=f"{work_dir}/run-{i:04d}/params.json")

                    # check if K0 file exists and the move it
                    filename = f'k0/{self.params["CONTROL"]["M"]}_{self.params["CONTROL"]["P"]}_{REs[i]:06.2f}.dat'
                    if os.path.isfile(filename):
                        os.system(f"mv {filename} {work_dir}/run-{i:04d}/k0.dat")
                    else:
                        continue

                    # run the code
                    M = N-len([f for f in os.listdir(queue_dir) if f.endswith('.txt')])
                    mp_worker(['./film-ns', f'./{work_dir}/run-{i:04d}/', i, N, my_rank, M])

                    break
                except subprocess.CalledProcessError:
                    # failed to rm queue file means the ith job is taken
                    continue
            # if we ever finish the loop if means that all the jobs are taken
            # and we are done
            if not ongoing:
                break
        comm.Barrier()

    def multi_run2(self, REs, CAs, Ms, Ps, rom="wr", uid=1):
        """
        Runs the code over all of the supplied parameters. The lengths must
        either match or one must have a single value
        """
        # get length and convert to lists if required
        N = len(REs)

        # get comms details
        comm = MPI.COMM_WORLD
        my_rank = comm.Get_rank()
        p = comm.Get_size()
        M = 0  # number of runs this task has to perform
        for i in range(N):
            if i % p != my_rank:
                continue
            M = M + 1

        # set up "queue"
        work_dir = f"py-work-{rom}-{uid:04d}"
        output_dir = f"py-out-{rom}-{uid:04d}"
        queue_dir = f"{work_dir}/queue"
        if my_rank == 0:
            os.system(f"mkdir -p {output_dir}")
            os.system(f"mkdir -p {work_dir}")
            if not os.path.isdir(f"mkdir -p {queue_dir}"):
                os.system(f"mkdir -p {queue_dir}")
                for i in range(N):
                    os.system(f"touch {queue_dir}/{i:04d}.txt")
        comm.Barrier()

        while True:
            # find next task
            ongoing = 0
            for i in range(N):
                try:
                    subprocess.run(["rm", f"{work_dir}/queue/{i:04d}.txt"],
                                   stderr=subprocess.DEVNULL,
                                   stdout=subprocess.DEVNULL,
                                   check=True)
                    ongoing = 1

                    # set up directories, executable and parameter file
                    os.system(f"mkdir {work_dir}/run-{i:04d}")
                    os.system(f"mkdir {work_dir}/run-{i:04d}/out")
                    os.system(f"mkdir {work_dir}/run-{i:04d}/dump")
                    os.system(f"cp film-ns {work_dir}/run-{i:04d}/film-ns")
                    self.params["CONTROL"]["rom"] = rom
                    self.set_params(REs[i], CAs[i])
                    self.params["CONTROL"]["M"] = Ms[i]
                    self.params["CONTROL"]["P"] = Ps[i]
                    self.dump_params(file=f"{work_dir}/run-{i:04d}/params.json")

                    # check if (static) K0 file exists and the move it
                    if self.params["CONTROL"]["strategy"] == "static":
                        filename = f'k0/{self.params["CONTROL"]["M"]}_{self.params["CONTROL"]["P"]}_{REs[i]:06.2f}.dat'
                        if os.path.isfile(filename):
                            os.system(f"mv {filename} {work_dir}/run-{i:04d}/k0.dat")
                        else:
                            continue

                    # run the code
                    M = N-len([f for f in os.listdir(queue_dir) if f.endswith('.txt')])
                    mp_worker(['./film-ns', f'./{work_dir}/run-{i:04d}/', i, N, my_rank, M])

                    break
                except subprocess.CalledProcessError:
                    # failed to rm queue file means the ith job is taken
                    continue
            # if we ever finish the loop if means that all the jobs are taken
            # and we are done
            if not ongoing:
                break
        comm.Barrier()

    def multi_run3(self, params_list, uid=1):
        """
        Runs the code over all of the supplied parameters. The lengths must
        either match or one must have a single value
        """
        # get length and convert to lists if required
        N = len(params_list)

        # get comms details
        comm = MPI.COMM_WORLD
        my_rank = comm.Get_rank()
        p = comm.Get_size()
        M = 0  # number of runs this task has to perform
        for i in range(N):
            if i % p != my_rank:
                continue
            M = M + 1

        # set up "queue"
        work_dir = f"py-work-{uid:04d}"
        output_dir = f"py-out-{uid:04d}"
        queue_dir = f"{work_dir}/queue"
        if my_rank == 0:
            os.system(f"mkdir -p {output_dir}")
            os.system(f"mkdir -p {work_dir}")
            if not os.path.isdir(f"mkdir -p {queue_dir}"):
                os.system(f"mkdir -p {queue_dir}")
                for i in range(N):
                    # The queue consists of a dummy file and the relevant
                    # params file
                    os.system(f"touch {queue_dir}/{i:04d}.txt")
                    self.params = params_list[i]
                    self.dump_params(f"{queue_dir}/params-{i:04d}.json")
        comm.Barrier()

        while True:
            # find next task
            ongoing = 0
            for i in range(N):
                try:
                    subprocess.run(["rm", f"{work_dir}/queue/{i:04d}.txt"],
                                   stderr=subprocess.DEVNULL,
                                   stdout=subprocess.DEVNULL,
                                   check=True)
                    ongoing = 1

                    # read parameters from file
                    self.load_params(f"{work_dir}/queue/params-{i:04d}.json")

                    # set up directories, executable and parameter file
                    os.system(f"mkdir {work_dir}/run-{i:04d}")
                    os.system(f"mkdir {work_dir}/run-{i:04d}/out")
                    os.system(f"mkdir {work_dir}/run-{i:04d}/dump")
                    os.system(f"cp film-ns {work_dir}/run-{i:04d}/film-ns")
                    self.dump_params(file=f"{work_dir}/run-{i:04d}/params.json")

                    # check if (static) K0 file exists and the move it
                    if self.params["CONTROL"]["strategy"] == "static":
                        pass
                        # filename = f'k0/{self.params["CONTROL"]["M"]}_{self.params["CONTROL"]["P"]}_{REs[i]:06.2f}.dat'
                        # if os.path.isfile(filename):
                        #     os.system(f"mv {filename} {work_dir}/run-{i:04d}/k0.dat")
                        # else:
                        #     continue

                    # run the code
                    M = N-len([f for f in os.listdir(queue_dir) if f.endswith('.txt')])
                    mp_worker(['./film-ns', f'./{work_dir}/run-{i:04d}/', i, N, my_rank, M])

                    break
                except subprocess.CalledProcessError:
                    # failed to rm queue file means the ith job is taken
                    continue
            # if we ever finish the loop if means that all the jobs are taken
            # and we are done
            if not ongoing:
                break
        comm.Barrier()

    def multi_read(self, REs, CAs, rom="benney", uid=1):
        """
        Reads the output from multirun
        """
        # get length and convert to lists if required
        try:
            N = len(REs)
        except TypeError:
            N = 1
            REs = [REs]
        try:
            M = len(CAs)
        except TypeError:
            M = 1
            CAs = [CAs]

        # esure equal length
        if N != M:
            if N == 1:
                REs = 0*CAs + REs[0]
                N = M
            elif M == 1:
                CAs = 0*REs + CAs[0]
            else:
                return []

        comm = MPI.COMM_WORLD
        my_rank = comm.Get_rank()
        if my_rank != 0:
            return

        # collect in the data
        print("COLLECTING DATA", flush=True)
        work_dir = f"py-work-{rom}-{uid:04d}"
        out_dir = f"py-out-{rom}-{uid:04d}"
        rates = np.zeros((N, 7))

        Kheader = "x d"
        Kdata = np.zeros([2**self.params["SOLVER"]["level"], N+2])
        x = np.linspace(0, 30, num=2**8+1)
        x = x[:-1] + 0.5*(x[1]-x[0])
        Kdata[:, 0] = x  # x
        Kdata[:, 1] = np.exp((np.cos(2*np.pi*(x-9) / 30) - 1)/0.01)  # w = 0.1

        Fheader = "Re Ca Hbe Ube Hwr Uwr Hns Uns"
        Fdata = np.zeros([N, 8])
        for i in range(N):
            print(f"{i}: {REs[i]}", flush=True)
            self.set_params(REs[i], CAs[i])

            # gain matrix
            try:
                Kheader = f"{Kheader} R{REs[i]}C{CAs[i]}"
                K = np.loadtxt(f"{work_dir}/run-{i:04d}/out/K.dat")
                Kdata[:, i+2] = K[1, :]
            except:
                pass

            # final data
            Fdata[i, 0] = REs[i]
            Fdata[i, 1] = CAs[i]
            try:
                F = np.loadtxt(f"{work_dir}/run-{i:04d}/out/benney-0.dat")
                Fdata[i, 2] = F[-1, 1]
                Fdata[i, 3] = F[-1, 3]
            except:
                pass
            try:
                F = np.loadtxt(f"{work_dir}/run-{i:04d}/out/wr-0.dat")
                Fdata[i, 4] = F[-1, 1]
                Fdata[i, 5] = F[-1, 3]
            except:
                pass
            try:
                F = np.loadtxt(f"{work_dir}/run-{i:04d}/out/ns-0.dat")
                Fdata[i, 6] = F[-1, 1]
                Fdata[i, 7] = F[-1, 3]
            except:
                pass

            # rates
            rates[i, 0] = REs[i]
            rates[i, 1] = CAs[i]
            try:
                self.read_data(model="benney", root=f"{work_dir}/run-{i:04d}")
                rates[i, 2] = self.rates["analytic_be"]
                rates[i, 3] = self.rates["analytic_wr"]
                rates[i, 4] = self.rates["benney"]
            except:
                pass

            try:
                self.read_data(model="wr", root=f"{work_dir}/run-{i:04d}")
                rates[i, 2] = self.rates["analytic_be"]
                rates[i, 3] = self.rates["analytic_wr"]
                rates[i, 5] = self.rates["wr"]
            except:
                pass

            try:
                self.read_data(model="ns", root=f"{work_dir}/run-{i:04d}")
                rates[i, 2] = self.rates["analytic_be"]
                rates[i, 3] = self.rates["analytic_wr"]
                rates[i, 6] = self.rates["ns"]
            except IndexError:
                pass

        subprocess.run(["rm", "-f", "out/Aloc.dat", "out/F.dat", "out/K.dat",
                        "out/J_be.dat", "out/J_wr.dat", "out/Oloc.dat",
                        "out/Phi_be.dat", "out/Phi_wr.dat", "out/Psi_be.dat",
                        "out/Psi_wr.dat"],
                       stderr=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL,
                       cwd=f'./{work_dir}/run-{i:04d}/')

        # save to a file
        M = self.params['CONTROL']['M']
        np.savetxt(f"{out_dir}/rates-{rom}-{M}.dat", rates)
        np.savetxt(f"{out_dir}/Ks-{rom}-{M}.dat", Kdata, header=Kheader, comments='')
        np.savetxt(f"{out_dir}/final-{rom}-{M}.dat", Fdata, header=Fheader, comments='')
        return rates

    def read_data(self, model="ns", root="."):
        """
        Reads in the data generated by run_model
        """
        # read the data
        try:
            self.data0[model] = np.loadtxt(f"{root}/out/{model}-0.dat")
        except FileNotFoundError:
            return
        self.data1[model] = []
        n = len(self.data0[model])
        for i in range(0, n):
            try:
                self.data1[model].append(np.loadtxt(f"{root}/out/{model}-1-{i:010d}.dat"))
            except OSError:
                n = i
                break

        # estimate damping rate from data
        i0 = 0
        i1 = n-1
        for i in range(0, n-1):
            if self.data0[model][i, 0] <= 0 and self.data0[model][i+1, 0] > 0:
                i0 = i
            if self.data0[model][i, 0] > 0 and self.data0[model][i, 1] < 1e-5:
                i1 = i
                if self.data0[model][i, 1] < 0:
                    i1 = i - 1
                break
        t = self.data0[model][i0:i1, 0]
        try:
            dh = np.log(self.data0[model][i0:i1, 1])
            m = len(t)
            X = np.ones([m, 2])
            X[0:m, 1] = t
            fit = np.linalg.lstsq(X, dh, rcond=None)
            self.rates[model] = fit[0][1]
        except OverflowError:
            self.rates[model] = 0

        # BENNEY
        # compute the linear damping rate
        J = np.loadtxt(f"{root}/out/A_be.dat")
        Psi = np.loadtxt(f"{root}/out/B_be.dat")
        K = np.loadtxt(f"{root}/out/K.dat")

        val, _ = np.linalg.eig(J - np.matmul(Psi, K))
        rate = np.max(np.real(val))
        self.rates["analytic_be"] = rate

        # create synthetic linear data
        if self.data0["analytic_be"] is None:
            self.data0["analytic_be"] = self.data0[model].copy()
        elif self.data0["analytic_be"][-1, 0] < self.data0[model][-1, 0]:
            self.data0["analytic_be"] = self.data0[model].copy()
        for i in range(0, len(self.data0["analytic_be"])):
            t = self.data0["analytic_be"][i, 0]
            if t < 0:
                self.data0["analytic_be"][i, 1] = 1
            else:
                try:
                    self.data0["analytic_be"][i, 1] = math.exp(rate*t)
                except OverflowError:
                    self.data0["analytic_be"][i, 1] = self.data0["analytic_be"][i-1, 1]

        # WR
        # compute the linear damping rate
        J = np.loadtxt(f"{root}/out/A_wr.dat")
        Psi = np.loadtxt(f"{root}/out/B_wr.dat")
        K = np.concatenate((K, np.zeros(K.shape)), axis=1)

        val, _ = np.linalg.eig(J - np.matmul(Psi, K))
        rate = np.partition(np.real(val).flatten(), -2)[-2]  # second largest
        self.rates["analytic_wr"] = rate

        # create synthetic linear data
        if self.data0["analytic_wr"] is None:
            self.data0["analytic_wr"] = self.data0[model].copy()
        elif self.data0["analytic_wr"][-1, 0] < self.data0[model][-1, 0]:
            self.data0["analytic_wr"] = self.data0[model].copy()
        for i in range(0, len(self.data0["analytic_wr"])):
            t = self.data0["analytic_wr"][i, 0]
            if t < 0:
                self.data0["analytic_wr"][i, 1] = 1
            else:
                try:
                    self.data0["analytic_wr"][i, 1] = math.exp(rate*t)
                except OverflowError:
                    self.data0["analytic_wr"][i, 1] = self.data0["analytic_wr"][i-1, 1]

    def plot_data(self, n=0):
        """
        Plots all recorded data
        """
        j = 0
        col = ["r", "y", "g", "b", "k"]
        for model in self.rates:
            if self.rates[model] is None:
                continue
            t = self.data0[model][:, 0]
            dh = np.exp(self.rates[model]*t)
            for i in range(0, len(t)):
                if t[i] >= 0:
                    break
                dh[i] = 1.0

            plt.plot(t, self.data0[model][:, 1], f"{col[j]}-", label=model)
            plt.plot(t, dh, f"{col[j]}--", label=f"{self.rates[model]:.2}")
            j = j + 1
        plt.gca().set_yscale('log')
        plt.ylim(top=5, bottom=1e-5)
        plt.xlabel("t")
        plt.ylabel("L2 deviation from target")
        plt.legend()
        plt.savefig(f"analysis-out/plot-{self.params['CONTROL']['rom']}-{n}-{self.Re:09.6f}.png")
        plt.clf()

def other_main():
    # get comms details
    comm = MPI.COMM_WORLD
    my_rank = comm.Get_rank()
    p = comm.Get_size()

    # Load params from the default file
    with open("params.json", 'r') as f:
        default_params = json.load(f)

    params_list = []
    n_dels = 40
    for i in range(n_dels):
        # for m in [5, 6, 7, 8, 9, 10, 11]:
        for m in [10]:
            DELs = np.linspace(0.0, 0.05 * m / 5.0, n_dels)

            params = copy.deepcopy(default_params)
            params["CONTROL"]["del"] = DELs[i]
            params["CONTROL"]["M"] = m
            params_list.append(params)

            params = copy.deepcopy(default_params)
            params["CONTROL"]["del"] = -DELs[i]
            params["CONTROL"]["M"] = m
            params_list.append(params)

    fm = FilmModel()
    fm.multi_run3(params_list)



def main(rom, n, uid, M):
    """
    runs across a range of Reynolds numbers
    """
    fm = FilmModel()

    Re0 = 1.0
    Re1 = 100.0
    Ca0 = 0.05
    Ca1 = 0.05

    Res = np.logspace(np.log10(Re0), np.log10(Re1), num=n)
    Ca = 0.05
    ms = [M]

    NN = len(Res)*len(ms)*4

    REs = np.zeros([NN])
    CAs = np.zeros([NN])
    Ms = np.zeros([NN], dtype=int)
    Ps = np.zeros([NN], dtype=int)

    k = 0
    for Re in Res:
        for m in ms:
            # for P in [m-1, m, m+1, m+2]:
            for P in [m]:
                REs[k] = Re
                CAs[k] = Ca
                Ms[k] = m
                Ps[k] = P
                k += 1

    fm.params["CONTROL"]["start"] = 50
    fm.params["DOMAIN"]["tmax"] = 400 + fm.params["CONTROL"]["start"]
    fm.params["SOLVER"]["dtout"] = 0.5
    fm.params["CONTROL"]["rom"] = rom
    fm.params["CONTROL"]["strategy"] = "lqr"  # TODO: rerun with dynamic

    # print(REs)
    # print(CAs)
    # print(Ms)
    # print(Ps)

    fm.multi_run2(REs, CAs, Ms, Ps, rom="wr", uid=uid)


if __name__ == "__main__":
    other_main()
    exit(0)

    # process args
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', action='store', type=int, default=10,
                        help='number of samples (default 10)')
    parser.add_argument('-u', '--uid', action='store', type=int,
                        help='unique ID for running multiple jobs (default M)')
    parser.add_argument('-m', action='store', type=int, default=5,
                        help='number of controls (default 5)')
    parser.add_argument('-w', '--wr', action='store_true',
                        help='only WR')
    parser.add_argument('-b', '--benney', action='store_true',
                        help='only Benney')
    args = parser.parse_args()

    if args.uid is None:
        args.uid = args.m

    # run requested code
    if args.wr:
        main("wr", args.n, args.uid, args.m)
    elif args.benney:
        main("benney", args.n, args.uid, args.m)
    else:
        main("wr", args.n, args.uid, args.m)
        main("benney", args.n, args.uid, args.m)
