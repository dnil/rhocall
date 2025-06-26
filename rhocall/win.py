class Win:
    """A general sequence window with one linear sum property.
    Outputs bed style, with a running average over the window.
    In this case, represents a run of autozygosity, with an average quality."""

    chr = ""
    start = 0
    end = 0
    sum = 0
    count = 0
    print_me = False

    def reset(self, chr="", start=0, end=0, sum=0, count=0, print_me=False):
        # current window object..
        self.chr = chr
        self.start = start
        self.end = end
        self.sum = sum
        self.count = count
        self.print_me = print_me

    def __init__(self, chr="", start=0, end=0, sum=0, count=0, print_me=False):
        self.reset(chr, start, end, sum, count, print_me)

    def extend(self, new_end, delta):
        self.end = new_end
        self.sum = self.sum + delta
        self.count = self.count + 1

    def dump_bed_header(self, output):
        output.write("#chr\tstart\tend\tAZ\tqual\n")

    def dump(self, output):
        if self.print_me:
            avg = float(self.sum / self.count)
            output.write("%s\t%d\t%d\t1\t%.2f\n" % (str(self.chr), self.start, self.end, avg))
