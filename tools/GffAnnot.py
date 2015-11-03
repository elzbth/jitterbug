class GffAnnot:


    def __init__(self, line):
                    (chrom, so_type, method, start, stop, score, strand, frame, tags) = line.split("\t")
                    self.chrom = chrom
                    self.type = so_type
                    self.method = method
                    self.start = int(start) - 1
                    self.stop = int(stop)
                    self.score = score
                    self.strand = strand
                    self.frame = frame
                    self.tags = {}
                    self.gff_line = line
                    self.length = (self.stop - self.start) + 1

                    tag_pairs = tags.split(";")
                    for pair in tag_pairs:
                                    #print pair
                                    pair.strip()
                                    if pair == "":
                                        continue
                                    (key, value) = pair.split("=")
                                    key = key.strip()
                                    value = value.strip()
                                    self.tags[key] = value

    def __eq__(self, other):
        if self.chrom != other.chrom:
            return False
        if self.type != other.type:
            return False
        if self.method != other.method:
            return False
        if self.start != other.start:
            return False
        if self.stop != other.stop:
            return False
        if self.score != other.score:
            return False
        if self.strand != other.strand:
            return False
        if self.frame != other.frame:
            return False
        for (key, value) in self.tags.iteritems():
            if key not in other.tags:
                return False
            elif other.tags[key] != value:
                return False
        return True

    def to_string(self):
        tag_string = ""
        for tag, value in self.tags.items():
            tag_string += "%s=%s;" % (tag, str(value)) 
        return "\t".join([self.chrom, self.type, self.method, str(self.start), str(self.stop), self.score, self.strand, self.frame, tag_string])

