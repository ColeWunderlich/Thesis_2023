#!/usr/bin/env -S python -u

import argparse
import math
import scipy.sparse as sparse
import scipy.io as sio
import numpy as np
import sys
import os
import multiprocessing as mp

from timeit import default_timer as timer
from datetime import datetime as time
from datetime import timedelta


class effective_gene_length(object):
    def __init__(self):
        self.gene_lengths = {}
        self.effective_lengths = {}
        self.read_length = 100

    def get_length(self, gene):
        if gene in self.gene_lengths:
            return self.gene_lengths[gene]
        else:
            return False

    def get_all_lengths(self):
        return self.gene_lengths

    def get_effective_length(self, gene):
        if gene in self.effective_lengths:
            return self.effective_lengths[gene]
        else:
            return False

    def get_all_effective_lengths(self):
        return self.effective_lengths

    def calc_effective_length(self, length):
        length = int(length) - self.read_length + 1
        if length > 0:
            return 1.0 / length
        else:
            return 0.0

    def is_wellformed(self, items):
        if items[0][0] == '#':
            return False
        if len(items) >= 2:
            return True
        else:
            return False

    def build(self, filename):
        self.source_file = filename
        try:
            f = open(self.source_file, "r")
        except:
            print("Error in loading gene lengths file: ", self.source_file)
            sys.exit(1)

        for line in f:
            items = line.strip().split('\t')
            if self.is_wellformed(items):
                self.gene_lengths[items[0]] = int(items[-1])
                self.effective_lengths[items[0]] = self.calc_effective_length(items[-1])
            else:
                print("Skipping line: ", line)

        f.close()


class effective_te_length(object):
    def __init__(self):
        self.te_lengths = {}
        self.effective_lengths = {}
        self.read_length = 100

    def get_length(self, te):
        if te in self.te_lengths:
            return self.te_lengths[te]
        else:
            return False

    def get_all_lengths(self):
        return self.te_lengths

    def get_effective_length(self, te):
        if te in self.effective_lengths:
            return self.effective_lengths[te]
        else:
            return False

    def get_all_effective_lengths(self):
        return self.effective_lengths

    def is_wellformed(self, items):
        if items[0][0] == '#':
            return False
        if len(items) != 9:
            return False
        if len(items[8].split('transcript_id')) != 2:
            return False
        else:
            return True

    def find_transcript_id(self, annot_attributes):
        transcript_id = False
        for attrib_pair in annot_attributes.split(';'):
            if len(attrib_pair.split('transcript_id')) > 1:
                transcript_id = attrib_pair.split('"')[1]
        return transcript_id

    def calc_effective_length(self, length):
        length = int(length) - self.read_length + 1
        if length > 0:
            return 1.0 / length
        else:
            return 0.0

    def build(self, filename):
        self.source_file = filename
        start_idx = 1
        try:
            f = open(filename, 'r')
        except:
            print("Error in loading TE gtf: ", self.source_file)
            sys.exit(1)

        for line in f:
            items = line.strip().split('\t')
            if not self.is_wellformed(items):
                print("Skipping line: ", line)
                continue
            else:
                annot_start = int(items[3])
                annot_stop = int(items[4])
                transcript_id = self.find_transcript_id(items[8])

                if not transcript_id:
                    print("Skipping line: ", line)
                    continue
                else:
                    if transcript_id not in self.te_lengths:
                        self.te_lengths[transcript_id] = 0
                    self.te_lengths[transcript_id] += (annot_stop - annot_start) + start_idx

        f.close()

        for te in self.te_lengths:
            self.effective_lengths[te] = self.calc_effective_length(self.te_lengths[te])


class dedup_header(object):
    def __init__(self):
        self.filename = ''
        self.cell_barcodes = []
        self.annotations = []
        self.tot_cell_barcodes = 0
        self.tot_annotations = 0

    def linecount_correct(self, dedup_lines):
        self.tot_annotations = int(dedup_lines[0])
        self.tot_cell_barcodes = int(dedup_lines[1])
        if len(dedup_lines) != self.tot_annotations + self.tot_cell_barcodes + 2:
            print("Error in header, line count mismatch. {} + {} + 2 != {}".format(self.tot_annotations,
                                                                                   self.tot_cell_barcodes,
                                                                                   len(dedup_lines)))
            sys.exit(1)
        else:
            return True

    def pare_annotations(self, annotations):
        self.pared_annotations = []
        for annot in annotations:
            pared_annot = annot.split('::')[0]
            self.pared_annotations.append(pared_annot)

    def build(self, filename):
        with open(filename) as f:
            dedup_lines = f.read().splitlines()

        if self.linecount_correct(dedup_lines):
            self.cell_barcodes = dedup_lines[-1 * self.tot_cell_barcodes:]
            self.annotations = dedup_lines[2:len(dedup_lines) - self.tot_cell_barcodes]
        self.pare_annotations(self.annotations)


class equivalence_classes(object):
    def __init__(self):
        self.header = dedup_header()
        self.unmoored_umis = {}
        self.anchored_counts = {}

    def process_equivalence_line(self, line):
        line = line.strip().split('\t')
        line = [int(x) for x in line]
        cell_barcode_index = line[0]
        tot_annotations = line[1]
        tot_umis = line[-1]
        cell_barcode = self.header.cell_barcodes[cell_barcode_index]
        equiv_annots = line[2:-1]

        return cell_barcode, tot_annotations, tot_umis, equiv_annots

    def grab_annotations(self, equiv_annots, tot_annotations, annot_mask):
        for annotation in equiv_annots:
            annot_mask[0,annotation] += 1.0 / tot_annotations
        annot_mask = annot_mask.tocsr()
        return annot_mask

    def build(self, filename, header_filename):
        self.header.build(header_filename)
        blank = sparse.lil_array(np.zeros(len(self.header.annotations)))

        with open(filename) as f:
            for line in f:
                cell_barcode, tot_annotations, tot_umis, equiv_annots = self.process_equivalence_line(line)

                if tot_annotations == 1:
                    if cell_barcode not in self.anchored_counts:
                        self.anchored_counts[cell_barcode] = blank.copy()
                    self.anchored_counts[cell_barcode][0,equiv_annots[0]] += tot_umis

                else:
                    if cell_barcode not in self.unmoored_umis:
                        self.unmoored_umis[cell_barcode] = []
                    annot_mask = self.grab_annotations(equiv_annots, tot_annotations, blank.copy())
                    for i in range(tot_umis):
                        self.unmoored_umis[cell_barcode].append(annot_mask)

        for cbc in self.anchored_counts:
            self.anchored_counts[cbc] = self.anchored_counts[cbc].tocsr()  # once per cell


class expectation_maximization(object):
    def __init__(self, maxiter, incl_init, incl_ongo, eff_lengths, numproc):
        self.min_step = 1.0
        self.max_step = 1.0
        self.inc_step = 4.0
        self.tolerance = 0.0001
        self.max_iterations = maxiter
        self.incl_init = incl_init
        self.incl_ongo = incl_ongo
        self.effective_lengths = eff_lengths
        self.num_processes = numproc

    def normalize(self, sparse_vect):
        tmp_sum = np.sum(sparse_vect.data)
        if tmp_sum > 0:
            norm = 1.0 / tmp_sum
        else:
            norm = 0.0
        return sparse_vect.multiply(norm)

    def compute_abundances(self, weight_vect, umi_masks, blank):
        count_vect = blank.copy()  # YYY
        for mask in umi_masks:
            tmp_cnt = mask.multiply(weight_vect)
            tmp_cnt = self.normalize(tmp_cnt)
            count_vect += tmp_cnt
        return count_vect

    def calc_step_weights(self, counts):
        weights = counts.multiply(self.effective_lengths)
        weights = self.normalize(weights)
        return weights

    def calc_diff_magnitude(self, new_weights, old_weights):
        difference = new_weights - old_weights
        magnitude = math.sqrt(np.sum((difference * difference).data))
        return difference, magnitude

    def calc_alphaS(self, step_magnitude, course_correction):
        alphaS = step_magnitude / course_correction
        alphaS = max(self.min_step, min(self.max_step, alphaS))
        return alphaS

    def check_min_max(self, alphaS):
        if alphaS == self.max_step:
            self.max_step = self.inc_step * self.max_step

        if 0 > self.min_step == alphaS:  # check coverage, when would this condition occur?
            self.min_step = self.inc_step * self.min_step

    def preprocess(self, unmoored_counts, count_vect):
        umi_masks = []
        for count in unmoored_counts:
            count_vect += count
            umi_masks.append(count != 0)

        return count_vect, umi_masks

    def run(self, unmoored_counts, initial_fixed_count, ongoing_fixed_count, blank):
        count_vect, umi_masks = self.preprocess(unmoored_counts, blank.copy())  # YYY
        if self.max_iterations == 0:
            return count_vect

        initial_weights = self.calc_step_weights(count_vect + initial_fixed_count)

        iteration_cntr = 0

        while iteration_cntr < self.max_iterations:
            iteration_cntr += 1
            initial_counts = self.compute_abundances(initial_weights, umi_masks, blank)
            first_step_weights = self.calc_step_weights(initial_counts + ongoing_fixed_count)

            first_step_diff, first_step_magnitude = self.calc_diff_magnitude(first_step_weights, initial_weights)
            if first_step_magnitude < self.tolerance:
                break

            first_step_counts = self.compute_abundances(first_step_weights, umi_masks, blank)
            second_step_weights = self.calc_step_weights(first_step_counts + ongoing_fixed_count)

            second_step_diff, second_step_magnitude = self.calc_diff_magnitude(second_step_weights, first_step_weights)
            if second_step_magnitude < self.tolerance:
                initial_weights = second_step_weights
                break

            velocity, velocity_magnitude = self.calc_diff_magnitude(second_step_diff, first_step_diff)
            if velocity_magnitude < self.tolerance:
                initial_weights = first_step_weights  # confirm with Ying's math, seems off
                break

            course_correction = math.sqrt(abs(np.sum((first_step_diff * velocity).data)))

            alphaS = self.calc_alphaS(first_step_magnitude, course_correction)

            new_weights = sparse.csr_array(np.fmax(0.0, initial_weights.todense()[0] + (2 * alphaS) * first_step_diff + (alphaS ** 2) * velocity))  # YYY

            if abs(alphaS - 1.0) > 0.01:
                try:
                    new_counts = self.compute_abundances(new_weights, umi_masks, blank)
                    new_weights = self.calc_step_weights(new_counts + ongoing_fixed_count)
                except:
                    print("Error in EMupdate")
                    raise

            self.check_min_max(alphaS)

            new_weights, initial_weights = initial_weights, new_weights

        if iteration_cntr > self.max_iterations:
            print("did not converge by " + str(self.max_iterations) + " iterations. \n")

        return self.compute_abundances(initial_weights, umi_masks, blank)

    def call_by_cell(self, anch_dict, unmo_dict, blank):
        final_counts = {}
        init_fixed_cnts = {}
        ongo_fixed_cnts = {}
        for cbc in anch_dict:
            if self.incl_init:
                init_fixed_cnt = anch_dict[cbc]
            else:
                init_fixed_cnt = blank.copy()

            if self.incl_ongo:
                ongo_fixed_cnt = anch_dict[cbc]
            else:
                ongo_fixed_cnt = blank.copy()

            init_fixed_cnts[cbc] = init_fixed_cnt
            ongo_fixed_cnts[cbc] = ongo_fixed_cnt

        # run will be done on each cbc for each thread in the pool
        mp.set_start_method('forkserver')
        pool = mp.Pool(self.num_processes)
        tbl_list = [(cbc, pool.apply_async(self.run, args=(unmo_dict[cbc], init_fixed_cnts[cbc], ongo_fixed_cnts[cbc], blank.copy()))) for cbc in anch_dict]

        full_results = [(cbc, n.get()) for (cbc, n) in tbl_list]

        pool.close()
        pool.join()

        for cbc, unmo_cnt in full_results:
            final_counts[cbc] = unmo_cnt + anch_dict[cbc]
            final_counts[cbc] = sparse.csr_array(final_counts[cbc])

        return final_counts


def read_options(parser):
    args = parser.parse_args()

    if not os.path.isfile(args.resfile):
        print("File error: ", args.resfile)
        sys.exit(1)

    if not os.path.isfile(args.hedfile):
        print("File error: ", args.hedfile)
        sys.exit(1)

    if not os.path.isfile(args.tefile):
        print("File error: ", args.tefile)
        sys.exit(1)

    if not os.path.isfile(args.genfile):
        print("File error: ", args.gtffile)
        sys.exit(1)

    args.argtext = "\n".join(("# ARGUMENTS LIST:",
                              "# name = %s" % args.prefix,
                              "# results file = %s" % args.resfile,
                              "# header file = %s" % args.hedfile,
                              "# gene file = %s " % args.genfile,
                              "# TE file = %s " % args.tefile,
                              "# processors = %s " % args.numproc,
                              "# initialize EM with anchored = %s" % args.yes_init,
                              "# weight every step of EM with anchored = %s" % args.yes_ongo,
                              "# number of iterations of EM = %s" % args.maxiter))
    return args


def prepare_parser():
    desc = "Wraps expectation-maximization UMI redistribution to deduplicated data."

    exmp = "Example: TEsolo --res dedup_results.txt --hed dedup_header.txt --gen gene_lengths.txt --tes TE_annotation.gtf "

    parser = argparse.ArgumentParser(prog='TEwrap', description=desc, epilog=exmp)

    parser.add_argument('--res', metavar='dedup-results-file', dest='resfile', type=str, required=True,
                        help='file with results from deduplication by equivalence class')
    parser.add_argument('--hed', metavar='dedup-header-file', dest='hedfile', type=str, required=True,
                        help='deduplication header file, gene annotations as gene_id::[info], te annotations as transcript_id::[info]')
    parser.add_argument('--gen', metavar='genic-lengths-file', dest='genfile', type=str, required=True,
                        help='file with gene annotation lengths by gene_id')
    parser.add_argument('--tes', metavar='TE-GTF-file', dest='tefile', type=str, required=True,
                        help='GTF file for transposable element annotations, lengths will be calc by transcript_id')
    parser.add_argument('--project', metavar='name', dest='prefix', default='TEwrapper_out',
                        help='Name of this project. DEFAULT: TEwrapper_out')
    parser.add_argument('--iter', metavar='iterations-of-EM', dest='maxiter', type=int, default=100,
                        help='Number of rounds of EM to do. For flat sum with no EM set to 0. DEFAULT: 100')
    parser.add_argument('--init', metavar='include-anchors-initially', dest='yes_init', type=bool, default=False,
                        help='True or False on whether anchored counts are used for initializing the EM. DEFAULT: False')
    parser.add_argument('--ongo', metavar='include-anchors-ongoing', dest='yes_ongo', type=bool, default=False,
                        help='True or False on whether anchored counts are used in every EM step. DEFAULT: False')
    parser.add_argument('--threads', metavar='number-of-processors', dest='numproc', type=int, default=50,
                        help="the number of processing threads for this run. DEFAULT: 50")

    return parser


def generate_annotation_lengths(processed_annots, gene_lengths, te_lengths):
    lengths = []
    for annot in processed_annots:
        lengths.append(te_lengths.get_effective_length(annot) or gene_lengths.get_effective_length(annot))

    lengths = sparse.csr_array(lengths)

    return lengths


def main():
    args = read_options(prepare_parser())
    print(args.argtext)
    print(time.now())
    start = timer()
    print("reading gene lengths ...")

    try:
        gene_lengths = effective_gene_length()
        gene_lengths.build(args.genfile)

    except:
        sys.stdout.write("Error extracting gene lengths from file.")
        sys.exit(1)

    print(time.now())
    stop = timer()
    print("Previous step duration:\t", timedelta(seconds=stop-start),"seconds")
    start=timer()
    print("extracting TE lengths ...")

    try:
        te_lengths = effective_te_length()
        te_lengths.build(args.tefile)

    except:
        print("Error extracting TE lengths from gtf file.")
        sys.exit(1)

    print(time.now())
    stop = timer()
    print("Previous step duration:\t", timedelta(seconds=stop-start),"seconds")
    start=timer()
    print("processing deduplication results ...")

    try:
        equivalences = equivalence_classes()
        equivalences.build(args.resfile, args.hedfile)

    except:
        print("Error processing deduplication results.")
        sys.exit(1)

    print(time.now())
    stop = timer()
    print("Previous step duration:\t", timedelta(seconds=stop-start),"seconds")
    start=timer()
    print("{} cells initially found".format(max(len(equivalences.unmoored_umis), len(equivalences.anchored_counts))))

    eff_lengths = generate_annotation_lengths(equivalences.header.pared_annotations, gene_lengths, te_lengths)
    blank = sparse.csr_array(np.zeros(max(eff_lengths.shape)))

    em_loops = expectation_maximization(args.maxiter, args.yes_init, args.yes_ongo, eff_lengths, args.numproc)
    final_cnt_dict = em_loops.call_by_cell(equivalences.anchored_counts, equivalences.unmoored_umis, blank)

    print(time.now())
    stop = timer()
    print("Previous step duration:\t", timedelta(seconds=stop-start),"seconds")
    start=timer()
    print("assembling and outputting experiment-wide data table ...")
    checker = False
    tblstack = False
    for cbc in equivalences.header.cell_barcodes:
        if cbc not in final_cnt_dict:
            continue
        if not checker:
            tblstack = final_cnt_dict[cbc]
            checker = True
        else:
            tblstack = sparse.vstack([tblstack, final_cnt_dict[cbc]])

    tblstack.data = np.rint(tblstack.data)
    tblstack = tblstack.astype(int)
    tblstack.eliminate_zeros()

    # output sparse matrix files
    with open("{}.annots".format(args.prefix), 'w') as fil_annot:
        fil_annot.write("\n".join(equivalences.header.annotations))

    with open("{}.cbcs".format(args.prefix), 'w') as fil_cbcs:
        fil_cbcs.write("\n".join(equivalences.header.cell_barcodes))

    sio.mmwrite("{}.mtx".format(args.prefix), tblstack)

    print(time.now())
    stop = timer()
    print("Previous step duration:\t", timedelta(seconds=stop-start),"seconds")
    start=timer()
    print("... run complete.")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt! \n")
        sys.exit(0)
