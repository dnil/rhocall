import argparse
import numpy
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
import sys


# compute the ratio of homozygous snps across the genome
def generate_bins(input_vcf, window, filter, mnv, minqual, rsid, minaf, aftag, maxaf, minsnv):
    bins = {}
    contigs = {}
    # collect the number of hetrozygous and homozygous snps for each bin
    for line in input_vcf:
        if line[0] == "#":
            if "##contig=<ID=" in line:
                contig_id = line.split("##contig=<ID=")[-1].split(",")[0]
                contig_len = int(line.split(",length=")[-1].split(">")[0])
                contigs[contig_id] = contig_len
            #                sys.stderr.write("DEBUG: contig id %s len %i\n"%(contig_id, contig_len))

            if "#CHROM" in line:
                for contig in contigs:
                    bins[contig] = numpy.zeros((int(math.ceil(contigs[contig] / float(window))), 2))
            continue
        if filter and not "\tPASS\t" in line:
            continue

        content = line.strip().split()
        # skip mnv calls
        if not mnv:
            if not len(content[3]) == 1 and not len(content[4]) == 1:
                continue

        # filter low quality variants
        ok = True
        try:
            quality = float(content[5])
            if quality < minqual:
                ok = False
        except:
            pass
        if not ok:
            continue

        # rrsid filter(if enabled)
        if rsid and not content[2].startswith("rs"):
            #            sys.stderr.write("DEBUG: drop non rsID variant.")
            continue

        # allele frequency filter
        if minaf:
            if ";{}=".format(aftag) in content[7]:
                try:
                    af = float(content[7].split(";{}=".format(aftag))[-1].split(";")[0])
                    if af < minaf or af > maxaf:
                        continue
                except:
                    print("Error parsing allele frequency.")
                    continue
            else:
                continue

        pos = int(math.floor(int(content[1]) / float(window)))
        if "1/1" in content[-1] or "0/0" in content[-1] or "1|1" in content[-1]:
            bins[content[0]][pos][1] += 1
        elif "./." in content[-1] or "./1" in content[-1]:
            pass
        else:
            bins[content[0]][pos][0] += 1

    # compute ratios
    for chromosome in bins:
        tmp_ratios = []
        for chrom_bins in bins[chromosome]:
            nr_snps = chrom_bins[:,0] + chrom_bins[:,1]
            tmp_ratios = numpy.where(nr_snps < minsnv, "NaN", chrom_array[:,1] / float(nr_snps))
        bins[chromosome] = numpy.array(tmp_ratios)
    return bins


# take the rhocall output and find all reported ROH
def extract_roh(rho):
    roh = {}

    for line in rho:
        if line[0] == "#":
            continue
        elif "CHROMOSOME_TOT" in line:
            continue

        if line[0] + line[1] == "RG":
            content = line.strip().split()
            if not content[2] in roh:
                # New chromosome
                roh[content[2]] = []
            # Mark start and end
            roh[content[2]].append([int(content[3]), int(content[4])])
    return roh


# create the plots, one for each chromosome, and print them to the assigned directory
def generate_plots(binned_zygosity, roh, window, pointsize, out_dir):
    for chromosome in binned_zygosity:
        if "GL" in chromosome:
            continue

        posvector = numpy.array(range(0, len(binned_zygosity[chromosome]))) * window / 1000.0

        plt.figure(1)
        plt.subplot(111)

        fraction = plt.scatter(
            posvector,
            binned_zygosity[chromosome],
            c="black",
            s=pointsize,
            alpha=0.5,
            marker="o",
            label="fraction",
        )
        median_fraction = numpy.median(
            binned_zygosity[chromosome][numpy.where(binned_zygosity[chromosome] > -1)]
        )
        (median,) = plt.plot(
            [0, max(posvector)],
            [median_fraction, median_fraction],
            linewidth=4,
            color="green",
            alpha=0.75,
            label="median fraction",
        )
        (ROH,) = plt.plot([0, 0], [-1, -1], linewidth=2, color="red", label="RHO")

        plt.legend([median, fraction, ROH], ["median fraction", "Fraction", "RHO"])
        if chromosome in roh:

            roh[chromosome] = numpy.array(roh[chromosome]) / 1000
            for r in roh[chromosome]:
                ycoord = [1.1, 1.1]
                plt.plot([r[0], r[1]], ycoord, linewidth=4, color="red")
        # plt.ylim(ymax = 3*median_coverage, ymin = 0)
        plt.title("RHO plot on chromosome {}".format(chromosome))

        plt.ylim(ymin=0)
        plt.ylim(ymax=1.2)
        plt.xlabel("Positions(Kb)")
        plt.ylabel("Fraction of homozygous snps")
        figure = plt.gcf()
        figure.set_size_inches(16, 10)
        plt.savefig("{}/{}.png".format(out_dir, chromosome), dpi=100)
        plt.close()


# create the plots, one for each chromosome, and print them to the assigned directory
def generate_wig(binned_zygosity, roh, window, outfile_basename):

    wigf = open(outfile_basename + ".wig", "w")
    wigf.write('track type=wiggle_0 description="Fraction of homozygous snps"\n')

    bedf = open(outfile_basename + ".bed", "w")
    bedf.write('track name=rhocall description="Regions of autozygosity"\n')
    for chromosome in binned_zygosity:
        if "GL" in chromosome:
            continue

        # IGV range 0-1, mid 0 works fine.
        wigf.write("fixedStep chrom=%s start=1 step=%i\n" % (chromosome, window))
        for z in binned_zygosity[chromosome]:
            wigf.write("{}\n".format(z))

        if chromosome in roh:
            for r in roh[chromosome]:
                bedf.write("{}\t{}\t{}\n".format(chromosome, r[0], r[1]))
