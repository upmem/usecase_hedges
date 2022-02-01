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

dpuIntra = True
singleShot = False
# code.load_dpus()
test_dpu_encoder = False
test_dpu_decoder = True

coderates = array([NaN, 0.75, 0.6, 0.5, 1./3., 0.25, 1./6.]
                  )  # table of coderates 1..6

# user-settable parameters for this test
coderatecode = 3  # test this coderate in coderaetes table above
npackets = 20  # number of packets (of 255 strands each) to generate and test
totstrandlen = 300  # total length of DNA strand

strandIDbytes = 2  # ID bytes each strand for packet and sequence number
strandrunoutbytes = 2  # confirming bytes end of each strand (see paper)
hlimit = 5000  # maximum size of decode heap, see pape
leftprimer = "TCGAAGTCAGCGTGTATTGTATG"
# for direct right appending (no revcomp)
rightprimer = "TAGTGAGTGCGATTAAGCGTGTT"

# this test generates substitution, deletion, and insertion errors
# sub,del,ins rates to simulate (as multiple of our observed values):
(srate, drate, irate) = 1.5 * array([0.0238, 0.0082, 0.0039])

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
print "real CR for hedges (inner) level ", 1/((totstrandlen * 2.0) / (bytesperstrand * 8))
print ''
# payload bytes per packet of 255 strands
messbytesperpacket = strandsperpacket * messbytesperstrand
# set code rate with left and right primers
print("setcoderates : leftprimer ", leftprimer, " rprimer ", rightprimer)
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


maxpacket = 5


def dnatomess_dpu(dnapacket, decoded_reference):
    # HEDGES decode strands of DNA (assumed ordered by packet and ID number) to a packet
    baddecodes = 0
    erasures = 0
    mpacket = zeros([strandsperpacket, bytesperstrand], dtype=uint8)
    # everything starts as an erasure
    epacket = ones([strandsperpacket, bytesperstrand], dtype=uint8)

    code.load_dpu_decoder()
    code.push_global_params_dpus()
    print dnapacket[0:maxpacket, :]
    (errcode, mess_, _, _, _, _, perf_cycles) = code.decode_DPU(
        dnapacket[0:maxpacket, :], 8*bytesperstrand, 1)
    code.stop_dpus()
    code.load_dpu_decoder()
    code.push_global_params_dpus()
    (errcode, mess_, _, _, _, _, perf_inst) = code.decode_DPU(
        dnapacket[0:maxpacket, :], 8*bytesperstrand, 0)
    #    dnapacket[i, :], 8*bytesperstrand)

    code.stop_dpus()

    print "perf_cycles ", perf_cycles
    print "perf_inst ", perf_inst
    if perf_cycles > 0:
        print "inst/cycles ", 0.1 * perf_inst/perf_cycles
    #  exit(0)
    # print mess_
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
        print decoded_reference[i, :]
        print mpacket[i, :]
        mse = sum(abs(mpacket[i, :] - decoded_reference[i, :])**2)
        print 'DECODER MSE [CPU/DPU] ON SINGLE STRAND ', mse
        raw_input("Press Enter to continue...")
    return (mpacket, epacket, baddecodes, erasures)


def dnatomess(dnapacket):
    # HEDGES decode strands of DNA (assumed ordered by packet and ID number) to a packet
    baddecodes = 0
    erasures = 0

    mpacket = zeros([strandsperpacket, bytesperstrand], dtype=uint8)
    # everything starts as an erasure
    epacket = ones([strandsperpacket, bytesperstrand], dtype=uint8)

    for i in range(strandsperpacket):
        # print(dnapacket[i, :])
        if i + 1 > maxpacket:
            break
        (errcode, mess, _, _, _, _) = code.decode(
            dnapacket[i, :], 8*bytesperstrand)
        if errcode > 0:
            baddecodes += 1
            erasures += max(0, messbytesperstrand-len(mess))
        lenmin = min(len(mess), bytesperstrand)
        mpacket[i, :lenmin] = mess[:lenmin]
        epacket[i, :lenmin] = 0
    return (mpacket, epacket, baddecodes, erasures)

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


def extractplaintext(cpacket):
    # extract plaintext from a corrected packet
    plaintext = zeros(strandsperpacketmessage*messbytesperstrand, dtype=uint8)
    for i in range(strandsperpacketmessage):
        plaintext[i*messbytesperstrand:(i+1)*messbytesperstrand] = (
            cpacket[i, strandIDbytes:strandIDbytes+messbytesperstrand])
    return plaintext

# function to create errors in a bag (or packet) of DNA strands


def createerrors(dnabag, srate, drate, irate):
    # for testing: create errors in a bag of strands
    (nrows, ncols) = dnabag.shape
    newbag = zeros([nrows, ncols], dtype=uint8)
    for i in range(nrows):
        dna = code.createerrors(dnabag[i, :], srate, drate, irate)
        lenmin = min(len(dna), ncols)
        newbag[i, :lenmin] = dna[:lenmin]
    return newbag


# DO THE TEST
print "for each packet, these statistics are shown in two groups:"
print "1.1 HEDGES decode failures, 1.2 HEDGES bytes thus declared as erasures"
print "1.3 R-S total errors detected in packet, 1.4 max errors detected in a single decode"
print "2.1 R-S reported as initially-uncorrected-but-recoverable total, 2.2 same, but max in single decode"
print "2.3 R-S total error codes; if zero, then R-S corrected all errors"
print "2.4 Actual number of byte errors when compared to known plaintext input"

print "load dpus"
# code.push_global_params_dpus()
print "push globals to dpus"
badpackets = 0
Totalbads = zeros(8, dtype=int)

print "npackets", npackets

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
        print "ENCODER MSE [CPU/DPU] ", mse

    # simulate errors in DNA synthesis and sequencing
    obspack = createerrors(dnapack, srate, drate, irate)

    # decode
    (dpacket, epacket, baddecodes, erasures) = dnatomess(
        obspack)  # decode the strands

    if test_dpu_decoder:
        # code.load_dpu_decoder()
        # code.push_global_params_dpus()
        (dpacket_dpu, epacket_dpu, baddecodes_dpu, erasures_dpu) = dnatomess_dpu(
            obspack,  dpacket)  # decode the strands
        # code.stop_dpus()
        print"host =>", dpacket
        print"dpu  =>", dpacket_dpu
        mse = sum(abs(dpacket_dpu - dpacket)**2)
        print '[HEDGES][DECODER][PACKET][NSTRANDS=', maxpacket, ']', mse
        exit(0)
    raw_input("Press Enter to continue...")

    (cpacket, tot_detect, tot_uncorrect, max_detect, max_uncorrect,
     toterrcodes) = correctmesspacket(dpacket, epacket)

    # check against ground truth
    messcheck = extractplaintext(cpacket)
    badbytes = count_nonzero(messplain-messcheck)

    # print summary line
    Totalbads += array([baddecodes, erasures, tot_detect, max_detect,
                       tot_uncorrect, max_uncorrect, toterrcodes, badbytes])
    print("packet %3d: (bad Decode : %3d, reasures : %3d, tot_detect :\
%3d, max_detect : %3d) ( tot_uncorrect : %3d, max_uncorrect : \
%3d, totercodes : %3d, badbytes : %3d)" % (ipacket,
                                           baddecodes, erasures, tot_detect, max_detect, tot_uncorrect,
                                           max_uncorrect, toterrcodes, badbytes)),

    print("packet OK" if badbytes == 0 else "packet NOT ok")
    exit(0) if singleShot else None
    continue
    if badbytes:
        badpackets += 1
print("all packets OK" if not badpackets else "some packets had errors!")
print("TOT: (%4d %4d %4d %4d) (%4d %4d %4d %4d)" % tuple(Totalbads))
