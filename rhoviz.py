import argparse
import numpy
import  os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math

version = "42"

#compute the ratio of homozygous snps across the genome
def generate_bins(args):
    bins={}
    contigs={}
    #collect the number of hetrozygous and homozygous snps for each bin
    for line in open(args.input_vcf):
        if line[0] == "#":
            if "##contig=<ID=" in line:
                contig_id=line.split("##contig=<ID=")[-1].split(",")[0]
                contig_len=int(line.split(",length=")[-1].split(",")[0])
                contigs[contig_id]=contig_len

            if "#CHROM" in line:
                for contig in contigs:
                    bins[contig]=numpy.zeros((int(math.ceil(contigs[contig]/float(args.window))),2) )
            continue
        if not args.nofilter and not "\tPASS\t" in line:
            continue

        content=line.strip().split()
        #skip mnv calls
        if not args.mnv:
            if not len(content[3]) == 1 and not len(content[4]) == 1:
                continue

        #filter low quality variants
        ok=True
        try:
            quality=float(content[5])
            if quality < args.minqual:
                ok=False
        except:
            pass
        if not ok:
            continue

        #rrsid filter(if enabled)
        if args.rsid and not content[2].startswith("rs"):
            continue

        #allele frequency filter
        if args.minaf:
            if ";{}=".format(args.aftag) in content[7]:
                try: 
                    af=float(content[7].split(";{}=".format(args.aftag))[-1].split(";")[0])
                    if af < args.minaf or af > args.maxaf:
                        continue
                except:
                    print "error"
                    continue
            else:
                continue


        pos=int(math.floor(int(content[1])/float(args.window)))
        if "1/1" in content[-1] or "0/0" in content[-1]:
                bins[content[0]][pos][1]+=1
        elif "./." in content[-1] or "./1" in content[-1]:
            pass
        else:
                bins[content[0]][pos][0]+=1

    #compute ratios
    for chromosome in bins:
        tmp_ratios=[]
        for window in bins[chromosome]:
            if sum(window) < args.minsnv:
                tmp_ratios.append(-1)
            else:
                tmp_ratios.append(window[1]/float(window[1]+window[0]))
        bins[chromosome]=numpy.array(tmp_ratios)
    return  bins

# take the rhocall output and find all reported ROH
def extract_roh(args):
    roh={}
    chrom_tot={}
    for line in open(args.rho):
        if line[0] == "#":
            continue
        elif "CHROMOSOME_TOT" in line:
            continue

        content=line.strip().split()
        if not content[0] in roh:
            roh[content[0]]=[]
        roh[content[0]].append([int(content[1]),int(content[2])])         

    return(roh,chrom_tot)

#create the plots, one for each chromosome, and print them to the assigned directory
def generate_plots(binned_zygosity,roh,chrom_tot,args):
    for chromosome in binned_zygosity:
        if "GL" in chromosome:
            continue
        
        posvector=numpy.array( range(0,len(binned_zygosity[chromosome])))*args.window/1000.0

        plt.figure(1)
        plt.subplot(111)

        fraction=plt.scatter(posvector, binned_zygosity[chromosome], c="black",s=args.pointsize, alpha=0.5,marker = 'o',label='fraction')
        median_fraction=numpy.median( binned_zygosity[chromosome][ numpy.where( binned_zygosity[chromosome] > -1  ) ])
        median, =plt.plot([0,max(posvector)],[median_fraction,median_fraction], linewidth=4, color='green',alpha=0.75,label="median fraction")
        ROH, =plt.plot([0,0],[-1,-1], linewidth=2, color='red',label="RHO")
        
        plt.legend([median,fraction,ROH], ["median fraction","Fraction","RHO"])
        if chromosome in roh:
            
            roh[chromosome]=numpy.array(roh[chromosome])/1000
            for r in roh[chromosome]:
                ycoord=[1.1,1.1]
                plt.plot([ r[0] ,r[1] ], ycoord, linewidth=4, color='red')
        #plt.ylim(ymax = 3*median_coverage, ymin = 0)
        plt.title('RHO plot on chromosome {}'.format(chromosome))
        
        plt.ylim(ymin=0)
        plt.ylim(ymax=1.2)
        plt.xlabel('Positions(Kb)')
        plt.ylabel('Fraction of homozygous snps')
        figure = plt.gcf()
        figure.set_size_inches(16, 10)
        plt.savefig("{}/{}.png".format(args.dir,chromosome),dpi=100)
        plt.close()


def main(args):
    os.system( "mkdir {}".format(args.dir) )
    binned_zygosity=generate_bins(args)
    roh,chrom_tot=extract_roh(args)
    generate_plots(binned_zygosity,roh,chrom_tot,args)


## Argument parsing
parser = argparse.ArgumentParser(description='Call runs of autozygosity.')
parser.add_argument('--input_vcf','-i', type=str,required=True, help='Input (sorted) vcf file')
parser.add_argument('--dir','-d', type=str,required=True, help='output directory, the files will be named dir/chr.png, one picture is printed per chromosome')
parser.add_argument('--rho','-r', type=str,required=True, help='Input RHO file produced from rhocall')
parser.add_argument('--window','-w', type=int,default=10000, help='window size(bases)')
parser.add_argument('--minsnv','-m', type=int,default=2, help='minimum number of snvs for each plotted bin')
parser.add_argument('--maxsnv','-M', type=int,default=20, help='maximum number of snvs for each plotted bin')
parser.add_argument('--minaf',type=float,default=0.1, help='minimum allele frequency(this variable must be set to 0 if the allele frequency is not annotated)')
parser.add_argument('--maxaf',type=float,default=0.9, help='maximum allele frequency')
parser.add_argument('--mnv', action='store_true', help='include MNV')
parser.add_argument('--aftag', type=str,default="1000GAF", help='the allele frequency tag')
parser.add_argument('--minqual','-q', type=int,default=100, help='do not add snvs to a bin if there quality is less than this value')
parser.add_argument('--pointsize','-p', type=int,default=8, help='Size of the points (pixels)')
parser.add_argument('--rsid','-s', action='store_true', help='Skip variants not containing an rsid')
parser.add_argument('--nofilter','-n', action='store_true', help='include variants, even if they are not labeled PASS')
args = parser.parse_args()

main(args)
