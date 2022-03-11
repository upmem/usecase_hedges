# HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints
# William H. Press, John A. Hawkins, Stephen Knox Jones Jr, Jeffrey M. Schaub, Ilya J. Finkelstein
# submitted to Proceedings of the National Academy of Sciences
#
# Demonstration driver and installation validation program
#
# This file demonstrates the use of the separately compiled modules NRpyDNAcode and NRpyRS.
# (Those modules are compiled in C++ using Numerical Recipes.  See separate documentation.)
#
# We encode a specified number of packets from known plaintext, create DNA errors, then
# decode the DNA and compare.

import numpy
from numpy import *
import NRpyDNAcode as code
import NRpyRS as RS
import sys

import csv


test_dpu_encoder = False
test_dpu_decoder = True
test_dpu_statistics = False


coderates = array([NaN, 0.75, 0.6, 0.5, 1./3., 0.25, 1./6.]
                  )  # table of coderates 1..6

# user-settable parameters for this test
coderatecode = 5  # test this coderate in coderaetes table above
npackets = 1  # number of packets (of 255 strands each) to generate and test)
totstrandlen = 300  # total length of DNA strand

strandIDbytes = 2  # ID bytes each strand for packet and sequence number
strandrunoutbytes = 2  # confirming bytes end of each strand (see paper)

dpu_fake_packet_mul_factor = 2
maxpacket = 1000
dpu_profiling = 1
hlimit = 92500  # maximum size of decode heap, see pape
leftprimer = "TCGAAGTCAGCGTGTATTGTATG"
# for direct right appending (no revcomp)
rightprimer = "TAGTGAGTGCGATTAAGCGTGTT"

# this test generates substitution, deletion, and insertion errors
# sub,del,ins rates to simulate (as multiple of our observed values):
(srate, drate, irate) = 1.5 * array([0.0238, 0.0082, 0.0059])
#(srate, drate, irate) = 1.5 * array([0.0207, 0.0063, 0.0039])

# set parameters for DNA constrants (normally not changed, except for no constraint)
max_hpoly_run = 4  # max homopolymer length allowed (0 for no constraint)
GC_window = 12  # window for GC count (0 for no constraint)
max_GC = 8  # max GC allowed in window (0 for no constraint)
min_GC = GC_window - max_GC

# not normally user settable because assumed by Reed-Solomon outer code:
strandsperpacket = 255  # for RS(255,32)
strandsperpacketcheck = 32  # for RS(255,32)

# compute some derived parameters and set parameters in NRpyDNAcode module
leftlen = len(leftprimer)
rightlen = len(rightprimer)
strandlen = totstrandlen - leftlen - rightlen
strandsperpacketmessage = strandsperpacket - strandsperpacketcheck
(NSALT, MAXSEQ, NSTAK, HLIMIT) = code.getparams()  # get settable code parameters
code.setparams(8*strandIDbytes, MAXSEQ, NSTAK,
               hlimit)  # change NSALT and HLIMIT
bytesperstrand = int(strandlen*coderates[coderatecode]/4.)
messbytesperstrand = bytesperstrand - strandIDbytes - \
    strandrunoutbytes  # payload bytes per strand
# payload bytes per packet of 255 strands
messbytesperpacket = strandsperpacket * messbytesperstrand
# set code rate with left and right primers
# print("setcoderates : leftprimer ", leftprimer, " rprimer ", rightprimer)
code.setcoderate(coderatecode, leftprimer, rightprimer)
# set DNA constraints (see paper)
code.setdnaconstraints(GC_window, max_GC, min_GC, max_hpoly_run)

# define a source of plaintext bytes, either random or The Wizard of Oz in Esperanto
UseWiz = True
if UseWiz:
    wizoffset = 0
    wizfile = "data/WizardOfOzInEsperanto.txt"
    with open(wizfile, 'r') as myfile:
        wiztext = myfile.read()
    wizbytes = array([c for c in wiztext]).view(uint8)
    wizlen = len(wizbytes)

    def getwiz(n):  # return next n chars from wiztext
        global wizoffset, wizlen
        if wizoffset + n > wizlen:
            wizoffset = 0
        bytes = wizbytes[wizoffset:wizoffset+n]
        wizoffset += n
        return bytes
else:
    def getwiz(n):
        return random.randint(0, high=256, size=n, dtype=uint8)

# print "test source of plaintext: (below should be same Esperanto text twice - weird characters OK)"
# print wiztext[55722:55870]

# wizoffset = 55722
# print "".join([chr(x) for x in getwiz(148)])
# print ""

# functions to create sequential packets from the plaintext source, and R-S protect them


def createmesspacket(packno):  # packno in range 0..255 with value 2 for strandIDbytes
    packet = zeros([strandsperpacket, bytesperstrand], dtype=uint8)
    plaintext = zeros(strandsperpacketmessage*messbytesperstrand, dtype=uint8)
    for i in range(strandsperpacket):
        packet[i, 0] = packno  # note assumes value 2 for strandIDbytes
        packet[i, 1] = i
        if i < strandsperpacketmessage:
            ptext = getwiz(messbytesperstrand)
            packet[i, strandIDbytes:strandIDbytes+messbytesperstrand] = ptext
            plaintext[i*messbytesperstrand:(i+1)*messbytesperstrand] = ptext
    return (packet, plaintext)


def protectmesspacket(packetin):  # fills in the RS check strands
    packet = packetin.copy()
    regin = zeros(strandsperpacket, dtype=uint8)
    for j in range(messbytesperstrand):
        for i in range(strandsperpacket):
            regin[i] = packet[i, ((j+i) % messbytesperstrand)+strandIDbytes]
        regout = RS.rsencode(regin)
        for i in range(strandsperpacket):
            packet[i, ((j+i) % messbytesperstrand)+strandIDbytes] = regout[i]
    return packet

# functions to encode a packet to DNA strands, and decode DNA strands to a packet


def messtodna(mpacket):
    # HEDGES encode a message packet into strands of DNA
    filler = array([0, 2, 1, 3, 0, 3, 2, 1, 2, 0, 3, 1, 3, 1, 2,
                   0, 2, 3, 1, 0, 3, 2, 1, 0, 1, 3], dtype=uint8)
    dpacket = zeros([strandsperpacket, totstrandlen], dtype=uint8)
    for i in range(strandsperpacket):
        # print(" in python , encode ", len(mpacket[i, :]))
        dna = code.encode(mpacket[i, :])
        if len(dna) < totstrandlen:  # need filler after message and before right primer
            dnaleft = dna[:-rightlen]
            dnaright = dna[-rightlen:]
            dna = concatenate(
                (dnaleft, filler[:totstrandlen-len(dna)], dnaright))
            # n.b. this can violate the output constraints (very slightly at end of strand)
        dpacket[i, :len(dna)] = dna
    return dpacket


def messtodna_DPU(mpacket):
    # HEDGES encode a message packet into strands of DNA
    filler = array([0, 2, 1, 3, 0, 3, 2, 1, 2, 0, 3, 1, 3, 1, 2,
                   0, 2, 3, 1, 0, 3, 2, 1, 0, 1, 3], dtype=uint8)
    dpacket = zeros([strandsperpacket, totstrandlen], dtype=uint8)
    dpu_dnapacket = code.messtodna_DPU(mpacket)
    for i in range(strandsperpacket):
        dna = dpu_dnapacket[i]
        if len(dna) < totstrandlen:  # need filler after message and before right primer
            dnaleft = dna[: -rightlen]
            dnaright = dna[-rightlen:]
            dna = concatenate(
                (dnaleft, filler[:totstrandlen-len(dna)], dnaright))
            # n.b. this can violate the output constraints (very slightly at end of strand)
        dpacket[i, :len(dna)] = dna
    return dpacket


def meanEff(inst, cyc):
    inst = 1. * inst
    cyc = 1. * cyc
    e = numpy.divide(inst * 1., cyc)
    em = numpy.mean(e)
    c = numpy.sum(cyc)
    cm = numpy.mean(cyc)
    i = numpy.sum(inst)
    ci = numpy.mean(inst)
    return c, i, cm, ci, e, em


def dnatomess_dpu(dnapacket, decoded_reference, cpu_time):
    # HEDGES decode strands of DNA (assumed ordered by packet and ID number) to a packet
    baddecodes = 0
    erasures = 0
    mpacket = zeros([strandsperpacket, bytesperstrand], dtype=uint8)
    # everything starts as an erasure
    epacket = ones([strandsperpacket, bytesperstrand], dtype=uint8)

    # artificially double the number of strands to compensate
    # tasklet equilibrium
    # print dnapacket.shape
    for i in range(dpu_fake_packet_mul_factor - 1):
        dnapacket = concatenate((dnapacket,  dnapacket))

    code.load_dpu_decoder()
    code.push_global_params_dpus()
    (errcode, mess_, _, _, _, _,
        _NR_TASKLETS,
        _clocks_per_sec,
        _host_time,
        _total_cycles,
        _hypcompute_cycles,
        _io_cycles,
        _push_cycles,
        _pop_cycles,
        _decode_cycles,
        _hashfunc_cycles,
        _penality_cycles,
        _hypcompute_first_section_cycles,
        _hypload_cycles,
        _nbyte_written_heap,
        _nbyte_loaded_heap,
        _nbyte_written,
        _nbyte_loaded,
        _nb_cycles_total) = code.decode_DPU(
        dnapacket[0:maxpacket, :], 8*bytesperstrand, 2 if dpu_profiling else 0)
    code.stop_dpus()
    if dpu_profiling:
        code.load_dpu_decoder()
        code.push_global_params_dpus()
        (errcode, mess_, _, _, _, _,
            NR_TASKLETS,
            clocks_per_sec,
            host_time,
            total_cycles,
            hypcompute_cycles,
            io_cycles,
            push_cycles,
            pop_cycles,
            decode_cycles,
            hashfunc_cycles,
            penality_cycles,
            hypcompute_first_section_cycles,
            hypload_cycles,
            nbyte_written_heap,
            nbyte_loaded_heap,
            nbyte_written,
            nbyte_loaded,
            nb_cycles_total) = code.decode_DPU(
            dnapacket[0:maxpacket, :], 8*bytesperstrand, 1)
        code.stop_dpus()

        totalc,  totali, totalcm, totalim,  totale, totalem = meanEff(
            _total_cycles, total_cycles)
        pushc, pushi, pushcm, pushim, pushe, pushem = meanEff(
            _push_cycles, push_cycles)
        popc, popi, popcm, popim, pope, popem = meanEff(
            _pop_cycles, pop_cycles)
        hypcomputec, hypecomputei, hypcomputecm, hypcomputeim,  hypcomputee, hypcomputeem = meanEff(
            _hypcompute_cycles, hypcompute_cycles)
        ioc, ioi, iocm, ioim, ioe, ioem = meanEff(_io_cycles, io_cycles)
        decodec, decodei, decodecm, decodeim,  decodee, decodeem = meanEff(
            _decode_cycles, decode_cycles)
        hashfuncc, hashfunci, hashfuncm, hashfuncim,  hashfunce, hashfuncem = meanEff(
            _hashfunc_cycles, hashfunc_cycles)
        penalityc, penalityi, penalitycm, penalityim, penalitye, penalityem = meanEff(
            _penality_cycles, penality_cycles)
        hypcompute_first_sectionc, hypcompute_first_sectioni, hypcompute_first_sectioncm,  hypcompute_first_sectionim, hypcompute_first_sectione, hypcompute_first_sectionem = meanEff(
            _hypcompute_first_section_cycles, hypcompute_first_section_cycles)
        hyploadc, hyploadi, hyploadcm, hyploadim, hyploade, hyploadem = meanEff(
            _hypload_cycles, hypload_cycles)

        dpu_time = 1.0 * totalcm / clocks_per_sec
        dpu_effective_time = 1.0 * nb_cycles_total / clocks_per_sec
        eff_total = _nb_cycles_total / nb_cycles_total
        eff_effective_total = totalem

        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Tasklet Balancing % ", dpu_time / host_time
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Host Time\n{:.2f} secs".format(host_time)
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Dpu Time \n{:.2f}".format(dpu_time)
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Dpu Effective Time\n{:.2f} secs".format(dpu_effective_time)
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Dpu MegaCycles\n{:.2f}".format(totalcm * 1e-6)
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Dpu MegaCycles 2\n{:.2f}".format(nb_cycles_total * 1e-6)
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Dpu MegaInst\n{:.2f}".format(totalim * 1e-6)
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Dpu Effective Pipeline Eff % \n{:.2f}".format(eff_effective_total)
        print "[HEDGES][DECODER][DPU][PERF][SUMMARY] Dpu Pipeline Eff % \n{:.2f}".format(eff_total)

        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) IO (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(iocm * 1e-6, ioc, 100 * iocm / totalcm)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) push (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(pushcm * 1e-6, 100 * pushcm / totalcm)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) pop (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(popc * 1e-6, 100 * popc / totalc)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) decode (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(decodec * 1e-6, 100 * decodec / totalc)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) hashfunc (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(hashfuncc * 1e-6, 100 * hashfuncc / totalc)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) penality (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(penalityc * 1e-6, 100 * penalityc / totalc)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) hypcompute_first_section (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(hypcompute_first_sectionc * 1e-6, 100 * hypcompute_first_sectionc / totalc)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) hypload (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(hyploadc * 1e-6, 100 * hyploadc / totalc)
        print "[HEDGES][DECODER][DPU][PERF][CYCLES] (func) hypcompute (MegaCycle,  %) \n{:.2f} {:>10.2f} ".format(hypcomputec * 1e-6, 100 * hypcomputec / totalc)

        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) IO (Eff %) \n{:.2f} ".format(ioem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) push (Eff %) \n{:.2f}  ".format(pushem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) pop (Eff %) \n{:.2f} ".format(popem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) decode (Eff %) \n{:.2f} ".format(decodeem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) hashfunc (Eff %) \n{:.2f} ".format(hashfuncem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) penality (Eff %) \n{:.2f} ".format(penalityem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) hypcompute_first_section (Eff %) \n{:.2f} ".format(hypcompute_first_sectionem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) hypload (Eff %) \n{:.2f} ".format(hyploadem)
        print "[HEDGES][DECODER][DPU][PERF][PIPELINE_EFF] (func) hypcompute (Eff %) \n{:.2f} ".format(hypcomputeem)

        nbyte_written_heap = 1e-6 * nbyte_written_heap
        nbyte_loaded_heap = 1e-6 * nbyte_loaded_heap
        nbyte_written = 1e-6 * nbyte_written
        nbyte_loaded = 1e-6 * nbyte_loaded
        w_global = nbyte_written_heap + nbyte_written
        r_global = nbyte_loaded_heap + nbyte_loaded

        print "[DPU][TASKLET XFER ABSOLUTE] W (Mbytes) (global) \n{:.2f}  ".format(w_global)
        print "[DPU][TASKLET XFER ABSOLUTE] R (Mbytes) (global) \n{:.2f}  ".format(r_global)
        print "[DPU][TASKLET XFER ABSOLUTE] W (Mbytes) (heap) \n{:.2f}  ".format(nbyte_written_heap)
        print "[DPU][TASKLET XFER ABSOLUTE] R (Mbytes) (heap)  \n{:.2f} ".format(nbyte_loaded_heap)
        print "[DPU][TASKLET XFER ABSOLUTE] W (Mbytes) (core func) \n{:.2f} ".format(nbyte_written)
        print "[DPU][TASKLET XFER ABSOLUTE] R (Mbytes) (core func) \n{:.2f} ".format(nbyte_loaded)

        print "[DPU][TASKLET XFER BANDWIDTH] W (Mbytes/sec) (global) \n{:.2f}  ".format(w_global/dpu_time)
        print "[DPU][TASKLET XFER BANDWIDTH] R (Mbytes/sec) (global) \n{:.2f}  ".format(r_global/dpu_time)
        print "[DPU][TASKLET XFER BANDWIDTH] W (Mbytes/sec) (heap) \n{:.2f}  ".format(nbyte_written_heap/dpu_time)
        print "[DPU][TASKLET XFER BANDWIDTH] R (Mbytes/sec) (heap)  \n{:.2f} ".format(nbyte_loaded_heap/dpu_time)
        print "[DPU][TASKLET XFER BANDWIDTH] W (Mbytes/sec) (core func) \n{:.2f} ".format(nbyte_written/dpu_time)
        print "[DPU][TASKLET XFER BANDWIDTH] R (Mbytes/sec) (core func) \n{:.2f} ".format(nbyte_loaded/dpu_time)

        dimm_server = 20.
        dpu_cpu_acc = 1.*cpu_time / host_time
        rank_cpu_acc = (64. * cpu_time) / host_time
        dimm_cpu_acc = (128. * cpu_time) / host_time
        server_cpu_acc = (dimm_server*128. * cpu_time) / host_time
        Wdimm = 43.
        Wcpu = 50.

        dpu_decoding_throughput = dnapacket.shape[0] / host_time

        print "[DPU][TOTAL STRANDSLEN] {:.3f}".format(totstrandlen)
        print "[DPU][DECODING THROUGPOUT] (dna_sequences/sec/DPU) {:.3f} ".format(dpu_decoding_throughput)
        print "[DPU/CPU]  Acceleration {:.3f} ".format(dpu_cpu_acc)
        print "[RANK/CPU] Acceleration {:.3f} ".format(rank_cpu_acc)
        print "[DIMM/CPU] Acceleration {:.3f} ".format(dimm_cpu_acc)
        print "[{} DIMM CLOUD /CPU] Acceleration {:.3f} ".format(dimm_server, server_cpu_acc)

        print "[RANK/CPU] Energy Efficiency Gain {:.3f} ".format(rank_cpu_acc * 2 * (Wcpu/Wdimm))
        print "[DIMM/CPU] Energy Efficiency Gain {:.3f} ".format(dimm_cpu_acc * (Wcpu/Wdimm))
        print "[{} DIMM CLOUD /CPU] Energy Efficiency Gain {:.3f} ".format(dimm_server, server_cpu_acc * (Wcpu/(Wdimm)))
        hdecoder = ['CR', 'hlimit',  'srate', 'drate', 'irate']
        ddecoder = [[coderates[coderatecode],   hlimit, srate, drate, irate]]
        hbw = ['', 'total', 'heap', 'corefunc']
        dbw = [
            ['dpu R (total)',  r_global,  nbyte_loaded_heap, nbyte_loaded],
            ['dpu W (total)', w_global, nbyte_written_heap, nbyte_written],
            ['dpu R (MB/s)',  r_global/dpu_time, nbyte_loaded_heap /
                dpu_time,  nbyte_loaded/dpu_time],
            ['dpu W (MB/s)',  w_global/dpu_time, nbyte_written_heap /
                dpu_time, nbyte_written/dpu_time]
        ]
        hperf = ['', 'total', 'io', 'heap push', 'heap pop', 'decode',
                 'heap_hasging', 'penality', 'hypload', 'hypcompute']
        dperf = [
            ['host (s)', host_time],
            ['dpu Tasklet Balancing %', dpu_time/host_time],
            ['dpu % total', 1, iocm/totalcm, pushcm/totalcm, popcm / totalcm, decodecm/totalcm, hashfuncm/totalcm,
                penalitycm/totalcm, hyploadcm/totalcm, hypcomputecm/totalcm],
            ['dpu cycles', totalcm, iocm, pushcm,  popcm, decodecm, hashfuncm,
                penalitycm,  hyploadcm, hypcomputecm],
            ['dpu inst', totalim, ioim, pushim, popim, decodeim, hashfuncim,
                penalityim,  hyploadim, hypcomputeim],
            ['dpu pipeline efficiency', eff_total, 'NC', pushem, popem, decodeem, hashfuncem,
                penalityem,  hyploadem, hypcomputeem],
        ]
        hgain = ['', 'acc', 'energy']
        dgain = [
            ['total strand len', totstrandlen],
            ['dpu/cpu', '{:.3f}'.format(dpu_cpu_acc), ],
            ['rank/cpu', '{:.3f}'.format(rank_cpu_acc),
             '{:.3f}'.format(rank_cpu_acc * 2 * (Wcpu/Wdimm))],
            ['dimm/cpu',   '{:.3f}'.format(dimm_cpu_acc),
             '{:.3f}'.format(dimm_cpu_acc * 1.0 * (Wcpu/Wdimm))],
            ['clouddimm/cpu',  '{:.3f}'.format(server_cpu_acc),
             '{:.3f}'.format(server_cpu_acc * (Wcpu/Wdimm))],
            ['dpu decoding throughput (seq/sec/DPU)', '{:.3f}'.format(
                dpu_decoding_throughput)],
        ]

        if test_dpu_statistics:
            # runs without perf, for statiscit mode
            code.load_dpu_decoder()
            code.push_global_params_dpus()
            (errcode, mess_, _, _, _, _,
                NR_TASKLETS,
                clocks_per_sec,
                _,
                _,
                hypcompute_cycles,
                _,
                push_cycles,
                pop_cycles,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                _) = code.decode_DPU(
                dnapacket[0:maxpacket, :], 8*bytesperstrand, 0)
            code.stop_dpus()
        print '-----------CSV--------------'
        writer = csv.writer(sys.stdout)
        writer.writerow(hdecoder)
        writer.writerows(ddecoder)
        writer.writerow(hbw)
        writer.writerows(dbw)
        writer.writerow(hperf)
        writer.writerows(dperf)
        writer.writerow(hgain)
        writer.writerows(dgain)
        print '----------------------------'

    for i in range(strandsperpacket):
        if i + 1 > maxpacket:
            break
        # print 'dnatomes packet', i
        mess = mess_[i, :]
        if errcode > 0:
            baddecodes += 1
            erasures += max(0, messbytesperstrand-len(mess))
        lenmin = min(len(mess), bytesperstrand)
        mpacket[i, :lenmin] = mess[:lenmin]
        epacket[i, :lenmin] = 0
    return (mpacket, epacket, baddecodes, erasures)


def dnatomess(dnapacket):
    # HEDGES decode strands of DNA (assumed ordered by packet and ID number) to a packet
    baddecodes = 0
    erasures = 0
    host_time = 0
    mpacket = zeros([strandsperpacket, bytesperstrand], dtype=uint8)
    # everything starts as an erasure
    epacket = ones([strandsperpacket, bytesperstrand], dtype=uint8)

    for i in range(strandsperpacket):
        if i + 1 > maxpacket:
            break
        (errcode, mess, _, _, _, _, host_time_) = code.decode(
            dnapacket[i, :], 8*bytesperstrand)
        if errcode > 0:
            baddecodes += 1
            erasures += max(0, messbytesperstrand-len(mess))
        lenmin = min(len(mess), bytesperstrand)
        mpacket[i, :lenmin] = mess[:lenmin]
        epacket[i, :lenmin] = 0
        host_time = host_time + host_time_
    code.print_and_reset_reg_max_heap()
    return (mpacket, epacket, baddecodes, erasures, host_time)

# functions to R-S correct a packet and extract its payload to an array of bytes


def correctmesspacket(packetin, epacket):
    # error correction of the outer RS code from a HEDGES decoded packet and erasure mask
    packet = packetin.copy()
    regin = zeros(strandsperpacket, dtype=uint8)
    erase = zeros(strandsperpacket, dtype=uint8)
    tot_detect = 0
    tot_uncorrect = 0
    max_detect = 0
    max_uncorrect = 0
    toterrcodes = 0
    for j in range(messbytesperstrand):
        for i in range(strandsperpacket):
            regin[i] = packet[i, ((j+i) % messbytesperstrand)+strandIDbytes]
            erase[i] = epacket[i, ((j+i) % messbytesperstrand)+strandIDbytes]
        locations = array(argwhere(erase), dtype=int32)
        (decoded, errs_detected, errs_corrected,
         errcode, ok) = RS.rsdecode(regin, locations)
        tot_detect += errs_detected
        tot_uncorrect += max(0, (errs_detected-errs_corrected))
        max_detect = max(max_detect, errs_detected)
        max_uncorrect = max(max_uncorrect, max(
            0, (errs_detected-errs_corrected)))
        toterrcodes += (0 if errcode == 0 else 1)
        for i in range(strandsperpacket):
            packet[i, ((j+i) % messbytesperstrand)+strandIDbytes] = decoded[i]
    return (packet, tot_detect, tot_uncorrect, max_detect, max_uncorrect, toterrcodes)

# extract plaintext from a corrected packet


def extractplaintext(cpacket):
    plaintext = zeros(strandsperpacketmessage*messbytesperstrand, dtype=uint8)
    for i in range(strandsperpacketmessage):
        plaintext[i*messbytesperstrand:(i+1)*messbytesperstrand] = (
            cpacket[i, strandIDbytes:strandIDbytes+messbytesperstrand])
    return plaintext

# function to create errors in a bag (or packet) of DNA strands
# for testing: create errors in a bag of strands


def createerrors(dnabag, srate, drate, irate):
    (nrows, ncols) = dnabag.shape
    newbag = zeros([nrows, ncols], dtype=uint8)
    for i in range(nrows):
        dna = code.createerrors(dnabag[i, :], srate, drate, irate)
        lenmin = min(len(dna), ncols)
        newbag[i, :lenmin] = dna[:lenmin]
    return newbag


# DO THE TEST
print "[HEDGES TEST] for each packet, these statistics are shown in two groups:"
print "[HEDGES TEST] 1.1 HEDGES decode failures, 1.2 HEDGES bytes thus declared as erasures"
print "[HEDGES TEST] 1.3 R-S total errors detected in packet, 1.4 max errors detected in a single decode"
print "[HEDGES TEST] 2.1 R-S reported as initially-uncorrected-but-recoverable total, 2.2 same, but max in single decode"
print "[HEDGES TEST] 2.3 R-S total error codes; if zero, then R-S corrected all errors"
print "[HEDGES TEST] 2.4 Actual number of byte errors when compared to known plaintext input"
print "[HEDGES TEST] CR index {} : {}".format(coderatecode, coderates[coderatecode])
print "[HEDGES TEST] real inner CR for hedges (inner) level ", 1/((totstrandlen * 2.0) / (bytesperstrand * 8))
print '[HEDGES TEST] test_dpu_encoder', test_dpu_encoder
print '[HEDGES TEST] test_dpu_decoder', test_dpu_decoder
print '[HEDGES TEST] test_dpu_statistics', test_dpu_statistics
print '[HEDGES TEST] npackets', npackets
print '[HEDGES TEST] hlimit', hlimit
print '[HEDGES TEST] dpu_profiling', dpu_profiling

badpackets = 0
badpackets_ = 0
Totalbads = zeros(8, dtype=int)
Totalbads_ = zeros(8, dtype=int)


for ipacket in range(npackets):

    # encode
    messpack, messplain = createmesspacket(
        ipacket)  # plaintext to message packet
    rspack = protectmesspacket(messpack)  # Reed-Solomon protect the packet

    # (DPU) encode to strands of DNA containing payload messplain

    dnapack = messtodna(rspack)
    if test_dpu_encoder:
        code.load_dpu_encoder()
        code.push_global_params_dpus()
        dnapack_dpu = messtodna_DPU(rspack)
        code.stop_dpus()
        # encode to strands of DNA containing payload messplain
        mse = sum(abs(dnapack - dnapack_dpu)**2)
        print "[HEDGES][ENCODER][MSE (CPU - DPU)]", mse

    # simulate errors in DNA synthesis and sequencing
    obspack = createerrors(dnapack, srate, drate, irate)

    # decode
    (dpacket, epacket, baddecodes, erasures, host_time) = dnatomess(
        obspack)  # decode the strands

    cpu_time = host_time * dpu_fake_packet_mul_factor
    print "[CPU] (cpu_decoder) time [s] ", cpu_time

    (cpacket, tot_detect, tot_uncorrect, max_detect, max_uncorrect,
     toterrcodes) = correctmesspacket(dpacket, epacket)

    # check against ground truth
    messcheck = extractplaintext(cpacket)
    badbytes = count_nonzero(messplain-messcheck)

    # print summary line
    Totalbads += array([baddecodes, erasures, tot_detect, max_detect,
                       tot_uncorrect, max_uncorrect, toterrcodes, badbytes])
    print("[CPU] packet %3d: (bad Decode : %3d, reasures : %3d, tot_detect :\
%3d, max_detect : %3d) ( tot_uncorrect : %3d, max_uncorrect : \
%3d, totercodes : %3d, badbytes : %3d)" % (ipacket,
                                           baddecodes, erasures, tot_detect, max_detect, tot_uncorrect,
                                           max_uncorrect, toterrcodes, badbytes)),

    print("packet OK" if badbytes == 0 else "packet NOT ok")

    if test_dpu_decoder:
        (dpacket_dpu, epacket_dpu, baddecodes_dpu, erasures_dpu) = dnatomess_dpu(
            obspack,  dpacket, cpu_time)  # decode the strands
        print"host res (real packet) =>", dpacket
        print"dpu  res (real packet) =>", dpacket_dpu

        mse = sum(abs(dpacket_dpu - dpacket)**2)
        print '[HEDGES][DECODER][PACKET][NSTRANDS (WITH FAKE FOR PERF)=', dpu_fake_packet_mul_factor * dnapack.shape[0], ']', mse
        print "[HEDGES][ENCODER][MSE (CPU - DPU)]", mse
        # assert(mse == 0)

    if test_dpu_decoder:
        (cpacket_, tot_detect_, tot_uncorrect_, max_detect_, max_uncorrect_,
         toterrcodes_) = correctmesspacket(dpacket_dpu, epacket_dpu)

        # check against ground truth
        messcheck_ = extractplaintext(cpacket_)
        badbytes_ = count_nonzero(messplain-messcheck_)

        # print summary line
        Totalbads_ += array([baddecodes_dpu, erasures_dpu, tot_detect_, max_detect_,
                             tot_uncorrect_, max_uncorrect_, toterrcodes_, badbytes_])
        print("[DPU] packet %3d: (bad Decode : %3d, reasures : %3d, tot_detect :\
%3d,     max_detect : %3d) ( tot_uncorrect : %3d, max_uncorrect : \
%3d,     totercodes : %3d, badbytes : %3d)" % (ipacket,
                                               baddecodes_dpu, erasures_dpu, tot_detect_, max_detect_, tot_uncorrect_,
                                               max_uncorrect_, toterrcodes_, badbytes_)),

        print("[DPU] packet OK" if badbytes_ == 0 else "packet NOT ok")

    if badbytes:
        badpackets += 1
    if test_dpu_decoder:
        if badbytes_:
            badpackets_ += 1


print("all packets OK" if not badpackets else "some packets had errors!")
print("TOT: (%4d %4d %4d %4d) (%4d %4d %4d %4d)" % tuple(Totalbads))

if test_dpu_decoder:
    print("[DPU] all packets OK" if not badpackets_ else "some packets had errors!")
    print("[DPU] TOT: (%4d %4d %4d %4d) (%4d %4d %4d %4d)" % tuple(Totalbads_))
