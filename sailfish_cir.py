# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme: version 0.1 alpha , maybe, 
#


import os
import functools
import gffutils
import subprocess
import getopt

import random


__doc__ = ''' Sailfish-cir ver 0.11a
--------------
Usage: python sailfish_cir.py [options]

Options:

    -g  path to genomic sequence fasta file
    -a  path to gene annotation file, ie, .gtf or .gff files
    -r  path to single-end raw sequencing reads file.
    -1  path to raw pair-end reads, mate 1
    -2  path to raw pair-end reads, mate 2
    -c  path to CIRI output file to specify circular RNA
    --bed  path to bed file which contains circular RNA.
    -o  output folder that contains the index built by sailfish and quantification results
    -k  k-mer size used by sailfish to built index. default is 21
	--libtype   format string describing the library type of your reads. default is "IU", [read more on libtype of Sailfish](http://sailfish.readthedocs.org/en/master/library_type.html)
    --mll mean library length, this option is to fix up the effective length.
    -h/--help	print this help message


    other:
    python sailfish_cir_alpha.py convert foo.output
    this will convert this ciri output file into foo.bed and foo.mapping , the latter file contains information of circular RNA-host gene mapping
'''






#    todo: under development:
#    --genemap file to specify the transcript-gene mapping
#    --isoform switch whether to reveal the structure
#

__author__ = 'zerodel'

GFFREAD_CMD = "gffread"
SAILFISH_CMD = "sailfish"

WINDOW_WIDTH = 3


def tmpp(str_tmp):
    print("---> %s" % str_tmp)


# custom error defined here
class GTFerr(Exception):
    pass


class GTFitem_err(GTFerr):
    '''base class of exception of seq-file item utility
    '''
    pass


class ThisShouldBeAFolder(Exception):
    pass


class AttributionIncomplete(GTFitem_err):
    '''
    this happens when some attribution segment are lost in gtf file entry.
    '''
    pass


class NoGeneIDInGTF(AttributionIncomplete):
    pass


class NoTranscriptIDInGTF(AttributionIncomplete):
    pass


class SAMerr(Exception):
    pass


class NoAttributeStringInGTF(GTFitem_err):
    pass


class NoSuchFile(IOError):
    pass


class WrongArguments(Exception):
    pass


class FileTypeError(Exception):
    pass


# some fundamental classes. eg, gtf entry in gtf file, or ciri entry in CIRI output.
class GTFitem(object):
    """
    read single line of .seqfile file , and construct items ,
    """
    sampleAttribute = r'gene_id transcript_id exon_number gene_biotype gene_name p_id protein_id transcript_name tss_id'

    is_ensemble = True

    def __init__(self, line_in_gtf=""):
        """
        Constructor
        """

        if not line_in_gtf:   # the "null" condition
            self._seqname = ""
            self._source = ""
            self._feature = ''
            self._start = -1
            self._end = -1
            self._score = '.'
            self._strand = "+"
            self._frame = "."
            self.init_null_attribute()

        else:   # have content in the given gtf-file line
            self._parse_line(line_in_gtf)

    def seqname(self):
        return self._seqname

    def set_seqname(self, seqname):
        self._seqname = seqname

    def starts(self):
        return self._start

    def ends(self):
        return self._end

    def set_start(self, start):
        self._start = start

    def set_end(self, end):
        self._end = end

    def get_gene_id(self):
        return self._attributes["gene_id"]

    def get_transcript_id(self):
        return self._attributes["transcript_id"]

    def set_gene_id(self, new_id):
        self._attributes["gene_id"] = new_id

    def set_transcript_id(self, new_transcript_id):
        self._attributes["transcript_id"] = new_transcript_id

    def get_strand(self):
        return self._strand

    def get_attribute(self):
        return self._attributes

    def set_source(self, source):
        self._source = source

    def get_source(self):
        return self._source

    def set_feature(self, feature):
        self._feature = feature

    def get_feature(self):
        return self._feature

    def set_strand(self, strand):
        self._strand = strand

    def set_frame(self, frame):
        self._frame = frame

    def get_frame(self):
        return self._frame

    def _parse_line(self, line_in_gtf):
        """ parse a line in seq-file file ,
        only gene id and transcript id will be extracted from attribute string
        """
        element_gtf = line_in_gtf.strip().split("\t")

        gtf_attribute_string = element_gtf.pop()
        try:
            self._check_attribute_string(gtf_attribute_string)
        except Exception as e:
            raise e

        self._seqname = element_gtf.pop(0)
        self._source = element_gtf.pop(0)
        self._feature = element_gtf.pop(0)
        self._start = int(element_gtf.pop(0))
        self._end = int(element_gtf.pop(0))
        self._score = element_gtf.pop(0)
        self._strand = element_gtf.pop(0)
        self._frame = element_gtf.pop(0)
        self._attributes = self.attribute2dict(gtf_attribute_string)

    @staticmethod
    def _check_attribute_string(gtf_attribute):
        if not gtf_attribute:
            # if nothing in attribute string
            raise NoAttributeStringInGTF

        if "gene_id" not in gtf_attribute:
            raise NoGeneIDInGTF

        if "transcript_id" not in gtf_attribute:
            raise NoTranscriptIDInGTF

    @staticmethod
    def attribute2dict(gtf_attribute):
        """
        extract information from the attribute string of gtf file.

        :param gtf_attribute:
        :return:
        """
        return dict([(item.split()[0], item.split()[-1].strip('"'))
                     for item in gtf_attribute.strip().split(";") if item])

    def __eq__(self, other_gtf_item):
        return self._seqname == other_gtf_item.seqname() \
            and self._start == other_gtf_item.starts() \
            and self._end == other_gtf_item.ends() \
            and self._strand == other_gtf_item.get_strand() \
            and self.get_gene_id() == other_gtf_item.get_gene_id() \
            and self.get_transcript_id() == other_gtf_item.get_transcript_id()

    def __len__(self):
        return self._end - self._start + 1

    def init_null_attribute(self):
        """
        two mandatory attributes : gene_id and transcript_id
        :return:
        """
        self._attributes = dict()
        self._attributes.setdefault("gene_id", "")
        self._attributes.setdefault("transcript_id", "")

    def __str__(self):
        """
        the first eight element are following GFF format, and with a description of GTF
        :return:
        """
        return "\t".join([self._seqname,
                          self._source,
                          self._feature,
                          str(self._start),
                          str(self._end),
                          self._score,
                          self._strand,
                          self._frame,
                          self._attr2str()])

    def whether_ensemble(self, yes_or_no):
        self.is_ensemble = yes_or_no

    def _attr2str(self):
        attr_rebuild = "; ".join(['%s "%s"' % (key, self._attributes.get(key))
                                  for key in self.sampleAttribute.strip().split()
                                  if key in self._attributes.keys()])

        if self.is_ensemble:
            return "%s;" % attr_rebuild

        return attr_rebuild


# CIRI data parser definition here
CIRI_OUTPUT_FILE_HEADER = r'circRNA_ID      chr     circRNA_start   circRNA_end     #junction_reads SM_MS_SMS       #non_junction_reads     junction_reads_ratio    circRNA_type    gene_id junction_reads_ID'


class CIRIEntry(object):

    def __init__(self, string_line_in_ciri_output_format=""):
        """ construct an empty ciri entry or from a string.
        :param string_line_in_ciri_output_format: optional, a single string line in CIRI output file, except file header
        """
        if string_line_in_ciri_output_format:
            self._parse_line(string_line_in_ciri_output_format)
        else:
            self.id = ""
            self.chr = ""
            self.start = ""
            self.end = ""
            self.circRNA_type = ""
            self.gene_id = ""
            self.junction_reads = []

    def _parse_line(self, string_ciri):
        """:param string_ciri: a CIRI output file formatted string, except file header
        :return :None , set up your CIRIEntry object"""
        elements = string_ciri.strip().split("\t")
        self.junction_reads = elements.pop().split(",")
        self.gene_id = elements.pop()
        self.circRNA_type = elements.pop()
        self.id = elements.pop(0)
        self.chr = elements.pop(0)
        self.start = elements.pop(0)
        self.end = elements.pop(0)

    def __str__(self):
        return 'id:%s\nchr:%s\nstart:%s\nend:%s\ntype:%s\ngene:%s\n' % (
            self.id, self.chr, self.start, self.end, self.circRNA_type, self.gene_id
        )

    def to_dot_bed_string(self, remove_chr=False):
        """transfer this object into a .bed file string
        :param remove_chr :  boolean, since chr1 in UCSC is just 1 in Ensembl,
            this option decide whether should "chr" be removed
        """
        if remove_chr:
            chromosome_id = self.chr[3:]
        else:
            chromosome_id = self.chr

        return "\t".join([chromosome_id, self.start, self.end, self.id]).strip()


# .fasta file operations
class FastaEntry(object):
    def __init__(self):
        self.reset()

    def set_id(self, some_given_id):
        self._id = some_given_id

    def get_id(self):
        return self._id

    def is_adapter_added(self):
        return self._has_addapter

    def is_effective_length_fixed(self):
        return self._has_effective_length_fixed

    def reset(self):
        self._id = ""
        self._seq_string = ""
        self._seq_items = []
        self._has_addapter = False
        self._has_effective_length_fixed = False

    def add_seq_part(self, line):
        self._seq_items.append(line.strip())

    def shrink(self):
        if len(self._seq_items) > 0 and self._seq_items[0]:
            seq_line_stripped = [short_line.strip() for short_line in self._seq_items]
            self._seq_items = []
            self._seq_string = "".join(seq_line_stripped)
        else:
            pass

    def add_adapter(self, kmer_len=20):
        self.shrink()
        # here , we think if kmer_len is bigger than the whole sequence is acceptable
        if not self._has_addapter:
            self._seq_string = "%s%s" % (self._seq_string[-kmer_len:], self._seq_string)
            self._has_addapter = True
        else:
            pass

    def pad_for_effective_length(self, length_needed):
        self.shrink()

        #nts = "".join([random.choice("ACGT") for i in range(length_needed)])
        nts = "".join(["N" for i in range(length_needed)])

        if self._seq_string and  not self._has_effective_length_fixed:
            self._seq_string = "%s%s" % (nts, self._seq_string)
            self._has_effective_length_fixed = True
        else:
            pass

    def __str__(self):
        self.shrink()
        return "%s\n%s" % (self._id.strip(), self._seq_string.strip())


def transform_fasta(fa, tmp_name, method_of_transform_fa_entry=str):
    with open(tmp_name, "w") as transformed:
        with open(fa) as fasta_file_reader:
            coolie_fa_line = FastaEntry()
            while True:
                current_line = fasta_file_reader.readline().strip()
                if current_line:

                    if current_line.startswith(">"):
                        if coolie_fa_line.get_id():
                            transformed.write(method_of_transform_fa_entry(coolie_fa_line) + "\n")
                        coolie_fa_line.reset()
                        coolie_fa_line.set_id(current_line)
                    else:
                        if coolie_fa_line.get_id():
                            coolie_fa_line.add_seq_part(current_line)

                else:  # jump out the while loop
                    if coolie_fa_line.get_id():
                        transformed.write(method_of_transform_fa_entry(coolie_fa_line))
                    break


def add_adapter_fa(adapter_length):
    def adapter_in_front(fa_entry):
        fa_entry.add_adapter(kmer_len=adapter_length)
        return str(fa_entry)
    return adapter_in_front


def pad_for_effective_length(length_needed):
    def make_up_this_line(fa_entry):
        fa_entry.pad_for_effective_length(length_needed)
        return str(fa_entry)
    return make_up_this_line


def do_convert_in_site(fa, your_method=str):
    '''
     this step will replace the original fa file , use it with caution!!!
    :param fa:
    :param your_method:
    :return:
    '''
    import shutil
    name_formal, extent_formal = os.path.splitext(fa)
    tmp_name = name_formal + "_tmp" + extent_formal
    transform_fasta(fa, tmp_name, method_of_transform_fa_entry=your_method)
    shutil.move(tmp_name, fa)


def format_your_fasta(fasta):
    # check whether your fasta is safe for adding adapter
    class NotSingleLineSequence(Exception):
        pass

    def check_your_file(fasta_file):
        with open(fasta_file) as read_it:
            while True:
                seq_line_this = read_it.readline().strip()
                if seq_line_this:
                    if seq_line_this.startswith(">"):
                        seq_1 = read_it.readline()
                        this_line_suppose_to_be_a_name_for_seq = read_it.readline()
                        if not this_line_suppose_to_be_a_name_for_seq.startswith(">"):
                            raise NotSingleLineSequence
                else:
                    break

    # main part
    try:
        check_your_file(fasta)

    except NotSingleLineSequence:
        # ok , transform it
        do_convert_in_site(fasta)

    else:
        # do nothing
        print("ok ")


# todo: you need a interface for bed file, how to assign host gene for bed line
# ----> you need two file , one is bed, one is mapping ..
def transform_ciri_to_bed(ciri_output_file, output_bed_file="", isoform_gene_mapping_file=""):
    abs_ciri_dir = os.path.abspath(ciri_output_file)
    main_part_ciri_path = os.path.splitext(abs_ciri_dir)[0]
    if not output_bed_file:
        output_bed_file = ".".join([main_part_ciri_path, "bed"])

    if not isoform_gene_mapping_file:
        isoform_gene_mapping_file = ".".join([main_part_ciri_path, "mapping"])

    with open(ciri_output_file) as ciri_file:
            ciri_file.readline()  # file head should be skipped

            with open(output_bed_file, "w") as exporter:
                with open(isoform_gene_mapping_file, "w") as mapping_file:
                    for line in ciri_file:
                        ciri_line_entry = CIRIEntry(line.strip())

                        new_bed_line = ciri_line_entry.to_dot_bed_string()
                        mapping_string = export_mapping_of_circular_isoform(ciri_line_entry)

                        if new_bed_line and mapping_string:
                            exporter.write(new_bed_line + "\n")
                            mapping_file.write(mapping_string + "\n")
                        else:
                            # todo : leave a warning here ?
                            pass



def export_mapping_of_circular_isoform(some_ciri_entry):
    if some_ciri_entry.id and some_ciri_entry.gene_id:
        return "\t".join([some_ciri_entry.id, some_ciri_entry.gene_id])
    else:
        return ""


class PredictedCircularRegion(object):
    def __init__(self, args_tuple, **kwargs):

        if args_tuple:
            self.predict_id, self.seqid, self.start, self.end = args_tuple
            self.start = int(self.start)
            self.end = int(self.end)
        elif kwargs:
            self.seqid = kwargs.get("chromosome")
            self.start = int(kwargs.get("start"))
            self.end = int(kwargs.get("end"))
            self.predict_id = kwargs.get("given_id")
        else:
            raise WrongArguments

    def is_flanking(self, gff_feature):
        return self.start <= int(gff_feature.start) and self.end >= int(gff_feature.end)

    def extract_flanked_linear_entries(self, gffutils_database):
        # extract all isoform part from some database of gtf file .
        # here we assume all exon has attribution of 'transcript_id'
        transcript_exon_dict = {}
        for linear_isoform in gffutils_database.region(seqid=self.seqid, start=self.start, end=self.end,
                                                       featuretype="transcript"):
            corresponding_circular_exons = [exon for exon in
                                            gffutils_database.children(linear_isoform.id, featuretype="exon",
                                                                       order_by="start",
                                                                       limit=(self.seqid, self.start, self.end),
                                                                       completely_within=True)]
            if corresponding_circular_exons:
                transcript_exon_dict.setdefault(linear_isoform.id, corresponding_circular_exons)

        return transcript_exon_dict

    def permutate_flanking_exons(self, db):
        # permutate exons to get all possible circRNA isoforms.
        # todo : to be implemented
        pass

    @staticmethod
    def generate_exon_for_circular_isoform(host_seqname, start, end, host_gene_id, host_tran_id, strand="+", frame="."):
        artifical_exon = GTFitem()
        artifical_exon.set_start(int(start))
        artifical_exon.set_end(int(end))
        artifical_exon.set_gene_id(host_gene_id)
        artifical_exon.set_transcript_id(host_tran_id)
        artifical_exon.set_seqname(host_seqname)
        artifical_exon.set_source("ciri")
        artifical_exon.set_feature("exon")
        artifical_exon.set_strand(strand)
        artifical_exon.set_frame(frame)
        return artifical_exon

    def arrange_exons_the_naive_way(self, db):
        # find all exons and put theme together , that's all, and here we assume that region id follow pattern like "id@gene"

        exons_raw = list(set([(exon.seqid, exon.source, exon.start, exon.end, exon.strand, exon.frame) for exon in
                                 db.region(seqid=self.seqid, start=int(self.start), end=int(self.end), featuretype="exon")]))

        exon_filtered = []   # start filter exon objects
        for exon in exons_raw:
            exon_seqid, exon_source, exon_start, exon_end, exon_strand, exon_frame = exon

            if exon_seqid == 'chrM':
                continue

            if not exon_source in ["processed_transcript", "protein_coding"]:
                continue

            if exon_start < self.start - WINDOW_WIDTH:
                exon_start = self.start

            if exon_end > self.end + WINDOW_WIDTH:
                exon_end = self.end

            exon_filtered.append((exon_seqid, exon_source, exon_start, exon_end, exon_strand, exon_frame))


        artifical_exons = []
        if len(self.predict_id.split("@")) == 2:
            transcript_id, gene_id = self.predict_id.split("@")
        else:
            transcript_id = self.predict_id
            gene_id = "n/a"

        for exon_locus in exon_filtered:
            exon_seqid, exon_source, exon_start, exon_end, exon_strand, exon_frame = exon_locus
            artifical_exons.append(self.generate_exon_for_circular_isoform(host_seqname=exon_seqid,
                                                                           start=exon_start,
                                                                           end=exon_end,
                                                                           host_gene_id=gene_id,
                                                                           host_tran_id=transcript_id,
                                                                           strand=exon_strand,
                                                                           frame=exon_frame
                                                                            )
                                   )

        return artifical_exons


    def mark_extracted_exons(self, dict_transcript_exon):
        # this funciton is after the extract_flanked.... function

        marked_exons = []
        for transcript_id in dict_transcript_exon.keys():
            for exon in dict_transcript_exon[transcript_id]:
                neo_exon = simplify_this_feature(exon, new_source="circRNA",
                                                 new_transcript_id="%s@%s" % (self.predict_id, transcript_id))
                marked_exons.append(neo_exon)

        return marked_exons


def simplify_this_feature(feature_from_gffutils_db, new_source="", new_transcript_id=""):
    artifical_exon = GTFitem(str(feature_from_gffutils_db))
    formal_gene_id = artifical_exon.get_gene_id()
    formal_trans_id = artifical_exon.get_transcript_id()
    artifical_exon.init_null_attribute()
    artifical_exon.set_gene_id(formal_gene_id)
    artifical_exon.set_transcript_id(formal_trans_id)

    if new_source:
        artifical_exon.set_source(new_source)
    if new_transcript_id:
        artifical_exon.set_transcript_id(new_transcript_id)

    return artifical_exon


def mark_permutate_exon(some_feature):
    #todo: need to give more information in id......
    source_ciri = "Permutation"
    id_ciri = "circ_%s_%s_%s_%s" % (some_feature.chrom, str(some_feature.start), str(some_feature.end), some_feature.id)
    neo_feature = simplify_this_feature(some_feature, new_source=source_ciri, new_transcript_id=id_ciri)
    return neo_feature


def parse_bed_line(line):
    parts = line.strip().split()
    if len(parts) < 4:
        raise FileTypeError

    chr_name, start, end, isoform_id = parts[:4]
    return isoform_id, chr_name, start, end


def parse_ciri_line(line):
    parts = line.strip().split("\t")
    if len(parts) < 4:
        raise FileTypeError

    isoform_id, chr_name, start, end = parts[:4]
    return isoform_id, chr_name, start, end


def parse_ciri_as_region(ciri_output):
    with open(ciri_output) as ciri_reader:
        ciri_reader.readline()
        for line in ciri_reader:
            yield PredictedCircularRegion(parse_ciri_line(line))


def parser_bed_as_region(bed_output_no_header):
    with open(bed_output_no_header) as read_bed:
        for line in read_bed:
            yield PredictedCircularRegion(parse_bed_line(line))


def get_gff_database(gtf_file):

    path_main, file_part = os.path.split(gtf_file)

    file_body_name, file_suffix = file_part.split(".")

    if "gtf" == file_suffix:
        db_file_path = os.path.join(path_main, ".".join([file_body_name, "db"]))
        if os.path.exists(db_file_path):
            db = gffutils.FeatureDB(db_file_path)
        else:
            # todo: here lies some trap
            #  due to difference between versions of .gtf file, binary database building process may be time exhausting
            db = gffutils.create_db(gtf_file, db_file_path)
    elif "db" == file_suffix:
        db = gffutils.FeatureDB(gtf_file)
    else:
        raise NameError
    return db


def put_gtf_file_same_dir_with_prediction_file(circular_prediction_file):
    output_gtf_path_name = os.path.join(os.path.split(circular_prediction_file)[0],
                                        os.path.split(circular_prediction_file)[-1].split(".")[0] + ".gtf")
    return output_gtf_path_name


def do_make_gtf_for_circular_prediction(gff_db, circular_candidate_regions, output_gtf_path_name="", is_isoform_structure_shown=True):
    """
    this function produce a .gtf for all circular RNA candidates

    -----
    gff_db: an gffutils database object.
    circular_candidate_regions: a list of PredictedCircularRegion object.
    output_gtf_path_name: a string specify the output gtf file .
    is_isoform_structure_shown: whether to show isoform structures , default is True
    """
    with open(output_gtf_path_name, "w") as gtf_items_writer:
            for region_circular in circular_candidate_regions:

                if is_isoform_structure_shown:
                    flanked_linear_entries = region_circular.extract_flanked_linear_entries(gff_db)
                    exons_marked_circular = region_circular.mark_extracted_exons(flanked_linear_entries)
                else:
                    exons_marked_circular = region_circular.arrange_exons_the_naive_way(db=gff_db)

                for simple_exon_marked in exons_marked_circular:
                    gtf_items_writer.write("%s\n" % str(simple_exon_marked))


def exec_this(list_of_args):
    try:

        subprocess.check_call(list_of_args, stdout=subprocess.PIPE)

    except ImportError:
        print("unable to import subprocess")
        sys.exit(-1)

    except OSError:
        print("may be some executive file is missing in %s" % str(list_of_args))
        sys.exit(-1)

    except subprocess.CalledProcessError:
        print("error in executing the commands")
        sys.exit(-1)

    else:
        pass


def parse_parameters(cmd_args, short_option_definition, long_option_definition):
    try:
        opts, args = getopt.gnu_getopt(cmd_args, short_option_definition, long_option_definition)
    except ImportError as e:
        print("unable to import python module 'getopt'")
        print(e)
        sys.exit(-1)

    except getopt.GetoptError as e:
        print("error when parsing command line arguments")
        print(e)
        sys.exit(-1)

    else:
        return opts, args


def build_cmd_gffread(gff_file, genomic_seqs, output_fasta, transcript_filter="CME"):
    gff_file = os.path.abspath(gff_file)
    genomic_seqs = os.path.abspath(genomic_seqs)
    output_fasta = os.path.abspath(output_fasta)

    cmds = [GFFREAD_CMD, gff_file, "-g", genomic_seqs]

    if transcript_filter:
        cmds.append("-%s" % transcript_filter)
    else:
        pass

    cmds.append("-w")
    cmds.append(output_fasta)

    return cmds


def build_cmd_sailfish_index(ref_transcripts, out_dir, kmer_len=None):
    ref_transcripts = os.path.abspath(ref_transcripts)
    out_dir = os.path.abspath(out_dir)

    cmds = [SAILFISH_CMD, "index", "-t", ref_transcripts, "-o", out_dir]
    if kmer_len:
        cmds.append("-k")
        cmds.append(str(kmer_len))

    return cmds


def build_cmd_sailfish_quant(index_dir, libtype, unmated_seq="", mate1="", mate2="", quant_dir=".", geneMap=""):
    unmated_seq = os.path.abspath(unmated_seq)
    mate1 = os.path.abspath(mate1)
    mate2 = os.path.abspath(mate2)
    quant_dir = os.path.abspath(quant_dir)

    def basic_quant_cmd(index_dir, libtype):
        quant_cmds = [SAILFISH_CMD, "quant", "-i", index_dir, "-l", libtype]
        return quant_cmds

    def add_args_pair_ends_read(mate1, mate2, quant_cmds):
        quant_cmds.extend(["-1", mate1, "-2", mate2])
        return quant_cmds

    def add_args_single_end_read(quant_cmds, unmated_seq):
        quant_cmds.extend(["-r", unmated_seq])
        return quant_cmds

    def add_args_quant_dir(quant_cmds, quant_dir):
        quant_cmds.extend(["-o", quant_dir])
        return quant_cmds

    def add_args_gene_map(quant_cmds, gene_map_gtf):
        quant_cmds.extend(["--geneMap", gene_map_gtf])
        return quant_cmds


    quantification_cmds = basic_quant_cmd(index_dir, libtype)

    if mate1 and mate2:
        quantification_cmds = add_args_pair_ends_read(mate1, mate2, quantification_cmds)
    elif len(unmated_seq) > 1:
        quantification_cmds = add_args_single_end_read(quantification_cmds, unmated_seq)
    else:
        raise WrongArguments

    quantification_cmds = add_args_quant_dir(quantification_cmds, quant_dir)

    if geneMap:
        quantification_cmds = add_args_gene_map(quantification_cmds, geneMap)

    return quantification_cmds


def do_extract_classic_linear_transcript(gff, fasta, output):
    cmds = build_cmd_gffread(gff_file=gff, genomic_seqs=fasta, output_fasta=output)
    exec_this(cmds)


def do_extract_circular_transcript(gff, fa, output):
    cmd = build_cmd_gffread(gff, fa, output, "")
    exec_this(cmd)


def do_add_addapt(fa, fa_name_after_decoration, adapter_length):
    do_convert_in_site(fa, add_adapter_fa(adapter_length))
    os.rename(fa, fa_name_after_decoration)

def do_make_up_for_effective_length(fa_name_before, length_for_effective_length):
    do_convert_in_site(fa_name_before, pad_for_effective_length(length_for_effective_length))

def do_combine_files(file_1, file_2, file_output):
    # cmds = "cat %s %s > %s" % (linear_fa, cirular_fa, whole_fa)
    # os.system(cmds)
    with open(file_output, "w") as output_lines:
        with open(file_1) as read1:
            for line in read1:
                output_lines.write(line)
        with open(file_2) as read2:
            for line in read2:
                output_lines.write(line)


def do_make_index_sailfish(ref_transcripts, index_folder, kmer_len=None):
    if not os.path.exists(index_folder):
        os.makedirs(index_folder)

    if not os.path.isdir(index_folder):
        raise ThisShouldBeAFolder

    cmds = build_cmd_sailfish_index(ref_transcripts, index_folder, kmer_len=kmer_len)
    exec_this(cmds)


def do_quant_sailfish_single_end(index_dir, libtype, single_end_seq_read="", quant_dir=".", geneMap=""):
    cmds = build_cmd_sailfish_quant(index_dir=index_dir, libtype=libtype,
                                    unmated_seq=single_end_seq_read, quant_dir=quant_dir, geneMap=geneMap)
    tmpp(" ".join(cmds))
    exec_this(cmds)


def do_quant_sailfish_pair_end(index_dir, libtype, mate1="", mate2="", quant_dir=".", geneMap=""):
    cmds = build_cmd_sailfish_quant(index_dir, libtype, mate1=mate1, mate2=mate2, quant_dir=quant_dir, geneMap=geneMap)
    tmpp(" ".join(cmds))
    exec_this(cmds)


class PipeLineQuantification():
    def __init__(self, list_console_cmd):
        self._customized_parameter = functools.partial(parse_parameters,
                                                       short_option_definition="c:g:a:r:1:2:o:k:h",
                                                       long_option_definition=["libtype=", "bed=", "help", "geneMap=", "isoform", "mll="])

        if list_console_cmd:
            opts, args = self._customized_parameter(cmd_args=list_console_cmd)
            self._assign_args(opts)

    def _assign_args(self, opts):
        setting_map_of_opts = dict(opts)

        if "-h" in setting_map_of_opts or "--help" in setting_map_of_opts:
            print(__doc__)
            sys.exit(0)

        if "-g" in setting_map_of_opts and setting_map_of_opts["-g"]:
            self._genomic_seq = setting_map_of_opts["-g"]
        else:
            print("genomic sequence file missing , exiting .....")
            sys.exit(-1)

        if "-a" in setting_map_of_opts and setting_map_of_opts["-a"]:
            self._annotation_file = setting_map_of_opts["-a"]
        else:
            print("annotation file  missing , exiting .....")
            sys.exit(-1)

        if "-o" in setting_map_of_opts and setting_map_of_opts["-o"]:
            # do output
            self._output_folder = setting_map_of_opts["-o"]
        else:
            print("quant output path missing, exiting ......")
            sys.exit(-1)

        self._is_reads_paired_end = False
        if "-r" in setting_map_of_opts and setting_map_of_opts["-r"]:
            # this is un-mated rate, ie, single end reads
            self._single_end_read_file_path = setting_map_of_opts["-r"]

        elif "-1" in setting_map_of_opts and "-2" in setting_map_of_opts \
                and setting_map_of_opts["-1"] and setting_map_of_opts["-2"]:
            # yes , this is paired end operation
            self._is_reads_paired_end = True
            self._paired_end_reads_1_file_path = setting_map_of_opts["-1"]
            self._paired_end_reads_2_file_path = setting_map_of_opts["-2"]
        else:   # neither single end nor paired end
            print("illegal specification of sequencing reads file, exiting .....")
            sys.exit(-1)

        self._kmer_len = 21
        if "-k" in setting_map_of_opts and setting_map_of_opts["-k"]:
            try:
                self._kmer_len = int(setting_map_of_opts["-k"])
            except ValueError as e:
                print(str(e))
                sys.exit(-1)

        self._index_libtype = "IU"
        if "--libtype" in setting_map_of_opts and setting_map_of_opts["--libtype"]:
            self._index_libtype = setting_map_of_opts["--libtype"]

        # give the circular prediction result.
        self._is_circ_isoform_provided = False

        if "-c" in setting_map_of_opts and setting_map_of_opts["-c"]:
            self._is_circ_isoform_provided = True
            self._ciri_output = setting_map_of_opts["-c"]

        elif "--bed" in setting_map_of_opts and setting_map_of_opts["--bed"]:
            self._bed_file = setting_map_of_opts["--bed"]
            self._is_circ_isoform_provided = True

        else:
            # no CIRI means you have to do it with only linear transcription
            pass

        if "--geneMap" in setting_map_of_opts and setting_map_of_opts["--geneMap"]:
            self._gene_mapping = setting_map_of_opts["--geneMap"]
        else:
            self._gene_mapping = ""

        self.reveal_isoform_structure = False
        if "--isoform" in setting_map_of_opts:
            self.reveal_isoform_structure = True

        self.mean_library_length = 150
        self.need_to_pad = False
        if "--mll" in setting_map_of_opts:
            self.mean_library_length = int(round(float(setting_map_of_opts["--mll"])))
            self.need_to_pad = True

    def process_the_pipe_line(self):
        # this the main part of perform the pipeline.
        if self._is_circ_isoform_provided:
                self._cicular_pipeline()
        else:
            tmpp("do the basic pipeline")
            self._basic_linear_pipeline()

    def _basic_linear_pipeline(self, linear_genomic_transcript="linear_transcript.fa"):
        linear_transcript_fasta_file_path = os.path.join(self._output_folder, linear_genomic_transcript)
        sailfish_index_folder_only_linear_transcript = os.path.join(self._output_folder, "index_only_linear")
        sailfish_quant_folder_only_linear_transcript = os.path.join(self._output_folder, "quant_only_linear")

        if not os.path.exists(linear_transcript_fasta_file_path):
            tmpp("extract the linear transcript")
            do_extract_classic_linear_transcript(self._annotation_file, self._genomic_seq, linear_transcript_fasta_file_path)
            do_convert_in_site(linear_transcript_fasta_file_path)

        tmpp("do make index now ")
        do_make_index_sailfish(linear_transcript_fasta_file_path, sailfish_index_folder_only_linear_transcript, self._kmer_len)

        tmpp("do quantification now ....")
        if self._is_reads_paired_end:
            tmpp("quant for pair end reads")
            do_quant_sailfish_pair_end(sailfish_index_folder_only_linear_transcript,
                                       libtype=self._index_libtype,
                                       mate1=self._paired_end_reads_1_file_path,
                                       mate2=self._paired_end_reads_2_file_path,
                                       quant_dir=sailfish_quant_folder_only_linear_transcript)
        else:
            tmpp("quant for single end reads")
            do_quant_sailfish_single_end(sailfish_index_folder_only_linear_transcript,
                                         libtype=self._index_libtype,
                                         single_end_seq_read=self._single_end_read_file_path,
                                         quant_dir=sailfish_quant_folder_only_linear_transcript)

    def _cicular_pipeline(self, linear_genomic_transcript="linear_transcript.fa"):
        linear_transcript_fasta_file_path = os.path.join(self._output_folder, linear_genomic_transcript)
        circular_only_annotation_exon_file_path = os.path.join(self._output_folder, "circular_only.gtf")
        circular_only_transcript_raw_fasta_file_path = os.path.join(self._output_folder, "circular_only_transcript_raw.fa")
        decorated_circular_transcript_fasta_file = os.path.join(self._output_folder, "circular_only_transcript_decorated.fa")
        whole_transcript_file_path = os.path.join(self._output_folder, "transcriptome_with_circular.fa")

        sailfish_index_folder_circular = os.path.join(self._output_folder, "index_circular")
        sailfish_quant_folder_circular = os.path.join(self._output_folder, "quant_circular")

        if not os.path.exists(linear_transcript_fasta_file_path):
            tmpp("extract linear transcript ..")
            do_extract_classic_linear_transcript(self._annotation_file, self._genomic_seq, linear_transcript_fasta_file_path)
            do_convert_in_site(linear_transcript_fasta_file_path)

        # prepare circular transcriptome fa file.

        tmpp("make a gtf for circular RNA")
        if not os.path.exists(circular_only_annotation_exon_file_path):
            # $where_we_make_gtf_for_circular_isoforms
            if self._bed_file:
                do_make_gtf_for_circular_prediction(get_gff_database(self._annotation_file),
                                                    parser_bed_as_region(self._bed_file),
                                                    circular_only_annotation_exon_file_path,
                                                    self.reveal_isoform_structure)
            elif self._is_circ_isoform_provided:
                #do_make_gtf_for_circ_exon_use_CIRI_output(self._annotation_file, self._ciri_output, circular_only_annotation_exon_file_path)
                do_make_gtf_for_circular_prediction(get_gff_database(self._annotation_file),
                                                    parse_ciri_as_region(self._ciri_output),
                                                    circular_only_annotation_exon_file_path,
                                                    self.reveal_isoform_structure)
            else:
                raise NoSuchFile("no circular prediction file....")
        else:
            tmpp("already has a annotation file here")

        tmpp("extract circular RNA transcript")
        if not os.path.exists(decorated_circular_transcript_fasta_file):
            do_extract_circular_transcript(circular_only_annotation_exon_file_path, self._genomic_seq, circular_only_transcript_raw_fasta_file_path)
            do_convert_in_site(circular_only_transcript_raw_fasta_file_path)
            do_add_addapt(circular_only_transcript_raw_fasta_file_path, decorated_circular_transcript_fasta_file, self._kmer_len -1)

            if self.need_to_pad and self.mean_library_length > self._kmer_len - 1:
                # todo : pad some random seq to complete the length
                do_make_up_for_effective_length(decorated_circular_transcript_fasta_file,
                                                self.mean_library_length - self._kmer_len + 1)

        else:
            tmpp("already has a dot-fasta sequence file for circular RNA")

        if not os.path.exists(whole_transcript_file_path):
            do_combine_files(linear_transcript_fasta_file_path, decorated_circular_transcript_fasta_file, whole_transcript_file_path)
        else:
            tmpp("combined sequence file already exists")

        tmpp("make index for sailfish")
        do_make_index_sailfish(ref_transcripts=whole_transcript_file_path, index_folder=sailfish_index_folder_circular, kmer_len=self._kmer_len)

        if self._is_reads_paired_end:
            tmpp("do pair end reads quantification")
            do_quant_sailfish_pair_end(sailfish_index_folder_circular,
                                       libtype=self._index_libtype,
                                       mate1=self._paired_end_reads_1_file_path,
                                       mate2=self._paired_end_reads_2_file_path,
                                       quant_dir=sailfish_quant_folder_circular, geneMap=self._gene_mapping)

        else:
            tmpp("do single end reads quantification")
            do_quant_sailfish_single_end(sailfish_index_folder_circular,
                                         libtype=self._index_libtype,
                                         single_end_seq_read=self._single_end_read_file_path,
                                         quant_dir=sailfish_quant_folder_circular, geneMap=self._gene_mapping)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print(__doc__)
    else:
        if sys.argv[1] == "convert":
            # here are the translation part.
            # todo : the convert part
            transform_ciri_to_bed(sys.argv[-1])

        else:
            work_on_it = PipeLineQuantification(sys.argv[1:])
            work_on_it.process_the_pipe_line()
