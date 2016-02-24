# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme: version 0.1 alpha , maybe, 
#


import os
import functools

__doc__ = ''' Sailfish-cir ver 0.1
--------------
Usage: python sailfish_cir.py [options]

Options:

    -g  path to genomic sequence fasta file
    -a  path to genomic annotation file, ie, .gtf or .gff files
    -r  path to your single-end sequencing reads file.
    -1  path to your pair-end reads , mate 1
    -2  path to your pair-end reads, mate 2
    -o  output folder that contains your quantification results
    -k  K-mer size used by sailfish to built index.
    --libtype   format string describing the library type of your reads.
    -h/--help   print this help message
    -c  path to CIRI output file to specify circular RNA, if provided, the -o result will contain estimation of circular transcript as well as linear ones.
'''

__author__ = 'zerodel'

GFFREAD_CMD = "gffread"
SAILFISH_CMD = "sailfish"


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


class No_gene_id_in_gtf(AttributionIncomplete):
    pass


class No_transcript_id_in_gtf(AttributionIncomplete):
    pass


class SAMerr(Exception):
    pass


class No_Attribute_Sting_in_gtf(GTFitem_err):
    pass


class NoSuchFile(IOError):
    pass


class WrongArguments(Exception):
    pass


# some fundamental classes. eg, gtf entry in gtf file, or ciri entry in CIRI output.
class GTFitem(object):
    """
    read single line of .seqfile file , and construct items ,
    """
    sampleAttribute = r'gene_id transcript_id exon_number gene_biotype gene_name p_id protein_id transcript_name tss_id'

    is_ensemble = True

    def __init__(self, line_in_gtf="", type_filter=""):
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
            self._strand = "."
            self._frame = "."
            self._init_null_attribute()

        else:   # have content in the given gtf-file line
            self._parse_line(line_in_gtf, type_filter="")

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

    def set_feature(self, feature):
        self._feature = feature

    def get_feature(self):
        return self._feature

    def set_strand(self, strand):
        self._strand = strand

    def _parse_line(self, line_in_gtf, type_filter=""):
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

    def _check_attribute_string(self,gtf_attribute):
        if not gtf_attribute:
            # if nothing in attribute string
            raise No_Attribute_Sting_in_gtf

        if "gene_id" not in gtf_attribute:
            raise No_gene_id_in_gtf

        if "transcript_id" not in gtf_attribute:
            raise No_transcript_id_in_gtf

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

    def _init_null_attribute(self):
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
        :param remove_chr :  boolean, since chr1 in UCSC is just 1 in Ensembl, this option decide whether should "chr" be removed
        """
        if remove_chr:
            chromesome_id = self.chr[3:]
        else:
            chromesome_id = self.chr

        return "\t".join([chromesome_id, self.start, self.end, self.id]).strip()

# .fasta file operations
class FastaEntry():
    def __init__(self):
        self.reset()

    def setid(self, some_given_id):
        self._id = some_given_id

    def getid(self):
        return self._id

    def is_addapter_added(self):
        return self._has_addapter

    def reset(self):
        self._id = ""
        self._seq_string = ""
        self._seq_items = []
        self._has_addapter = False

    def add_seq_part(self, line):
        self._seq_items.append(line.strip())

    def flush(self):
        if len(self._seq_items) > 0 and self._seq_items[0]:
            seq_line_stripped = [short_line.strip() for short_line in self._seq_items]
            self._seq_items = []
            self._seq_string = "".join(seq_line_stripped)
        else:
            pass

    def add_adapter(self, kmer_len=20):
        self.flush()
        # here , we think if kmer_len is bigger than the whole sequence is acceptable
        if not self._has_addapter:
            self._seq_string = "%s%s" % (self._seq_string[-kmer_len:], self._seq_string)
            self._has_addapter = True
        else:
            pass

    def __str__(self):
        self.flush()
        return "%s\n%s" % (self._id.strip(), self._seq_string.strip())


def transform_fasta(fa, tmp_name, method_of_transform_fa_entry=str):
    with open(tmp_name, "w") as transformed:
        with open(fa) as fasta_file_reader:
            coolie_fa_line = FastaEntry()
            while True:
                current_line = fasta_file_reader.readline().strip()
                if current_line:

                    if current_line.startswith(">"):
                        if coolie_fa_line.getid():
                            transformed.write(method_of_transform_fa_entry(coolie_fa_line) + "\n")
                        coolie_fa_line.reset()
                        coolie_fa_line.setid(current_line)
                    else:
                        if coolie_fa_line.getid():
                            coolie_fa_line.add_seq_part(current_line)

                else:  # jump out the while loop
                    if coolie_fa_line.getid():
                        transformed.write(method_of_transform_fa_entry(coolie_fa_line))
                    break


def add_adapter_fa(fa_entry):
    fa_entry.add_adapter()
    return str(fa_entry)


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
    # check whether your fasta is safe for adding addapter
    class NotSingleLineSequence(Exception):
        pass

    def check_your_file(fasta):
        with open(fasta) as readit:
            while True:
                seq_line_this = readit.readline().strip()
                if seq_line_this:
                    if seq_line_this.startswith(">"):
                        seq_1 = readit.readline()
                        this_line_suppose_to_be_a_name_for_seq = readit.readline()
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


# todo: you need a interface for bed file
# GTF file preparation here
def do_make_gtf_for_circ_exon(gtf_file, ciri_output, output_gtf_path_name=""):
    try:
        import gffutils
    except Exception as e:
        print("fail to import package : gffutils")
        raise e

    def generate_exon_for_circular_isoform(exon_locus, host_gene_id, host_seqname, host_tran_id):
        artifical_exon = GTFitem()
        artifical_exon.set_start(exon_locus[0])
        artifical_exon.set_end(exon_locus[-1])
        artifical_exon.set_gene_id(host_gene_id)
        artifical_exon.set_transcript_id(host_tran_id)
        artifical_exon.set_seqname(host_seqname)
        artifical_exon.set_source("ciri")
        artifical_exon.set_feature("exon")
        artifical_exon.set_strand("+")
        return artifical_exon

    path_main , file_part = os.path.split(gtf_file)

    file_body_name, file_suffix = file_part.split(".")

    if "gtf" == file_suffix:
        db_file_path = os.path.join(path_main, ".".join([file_body_name, "db"]))
        if os.path.exists(db_file_path):
            db = gffutils.FeatureDB(db_file_path)
        else:
            db = gffutils.create_db(gtf_file, db_file_path)
    elif "db" == file_suffix:
        db = gffutils.FeatureDB(gtf_file)
    else:
        raise NameError

    if not output_gtf_path_name:
        output_gtf_path_name = os.path.join(os.path.split(ciri_output)[0], os.path.split(ciri_output)[-1].split(".")[0] + ".gtf")

    print("output file is :", output_gtf_path_name)

    with open(ciri_output, "r") as read_ciri:
        read_ciri.readline()

        with open(output_gtf_path_name, "w") as exporter:
            for line in read_ciri:
                ciri_item = CIRIEntry(line.strip())
                host_gene_id = ciri_item.gene_id
                host_tran_id = ciri_item.id
                # here , hg19 file seqname has 'chr' before the chromesome number . so get rid of it .
                host_seqname  = ciri_item.chr[3:]

                if host_seqname and host_gene_id and host_tran_id and ciri_item.start and ciri_item.end and 'exon' == ciri_item.circRNA_type:
                    pass
                else:
                    continue

                exons = sorted(list(set(
                        [(exon.start, exon.end) for exon in
                         db.region(seqid= host_seqname, start=int(ciri_item.start), end=int(ciri_item.end), featuretype="exon")])))

                for exon_locus in exons:
                    artifical_exon = generate_exon_for_circular_isoform(exon_locus, host_gene_id, host_seqname,
                                                                        host_tran_id)

                    exporter.write(str(artifical_exon) + "\n")

    return output_gtf_path_name


def exec_this(list_of_args):
    try:
        import subprocess
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
        import getopt
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


def build_cmd_sailfish_quant(index_dir, libtype, unmated_seq="", mate1="", mate2="", quant_dir="."):
    unmated_seq = os.path.abspath(unmated_seq)
    mate1 = os.path.abspath(mate1)
    mate2 = os.path.abspath(mate2)
    quant_dir = os.path.abspath(quant_dir)

    def basic_quant_cmd(index_dir, libtype):
        quant_cmds = [SAILFISH_CMD, "quant", "-i", index_dir, "-l", libtype]
        return quant_cmds

    def add_args_pair_ends_read(mate1, mate2, quant_cmds):
        quant_cmds.append("-1")
        quant_cmds.append(mate1)
        quant_cmds.append("-2")
        quant_cmds.append(mate2)
        return quant_cmds

    def add_args_single_end_read(quant_cmds, unmated_seq):
        quant_cmds.append("-r")
        quant_cmds.append(unmated_seq)
        return quant_cmds

    def add_args_quant_dir(quant_cmds, quant_dir):
        quant_cmds.append("-o")
        quant_cmds.append(quant_dir)
        return quant_cmds

    quant_cmds = basic_quant_cmd(index_dir, libtype)

    if mate1 and mate2:
        quant_cmds = add_args_pair_ends_read(mate1, mate2, quant_cmds)
    elif len(unmated_seq) > 1:
        quant_cmds = add_args_single_end_read(quant_cmds, unmated_seq)
    else:
        raise WrongArguments

    quant_cmds = add_args_quant_dir(quant_cmds, quant_dir)

    return quant_cmds


def do_extract_classic_linear_transcript(gff, fasta, output):
    cmds = build_cmd_gffread(gff_file=gff, genomic_seqs=fasta, output_fasta=output)
    exec_this(cmds)


def do_extract_circular_transcript(gff, fa, output):
    cmd = build_cmd_gffread(gff, fa, output, "")
    exec_this(cmd)


def do_add_addapt(fa):
    do_convert_in_site(fa, add_adapter_fa)


def do_combine_circular_fa(linear_fa, cirular_fa, whole_fa):
    # cmds = "cat %s %s > %s" % (linear_fa, cirular_fa, whole_fa)
    # os.system(cmds)
    with open(whole_fa, "w") as output_lines:
        with open(linear_fa) as read1:
            for line in read1:
                output_lines.write(line)
        with open(cirular_fa) as read2:
            for line in read2:
                output_lines.write(line)


def do_make_index_sailfish(ref_transcripts, index_folder, kmer_len=None):
    if not os.path.exists(index_folder):
        os.makedirs(index_folder)

    if not os.path.isdir(index_folder):
        raise ThisShouldBeAFolder

    cmds = build_cmd_sailfish_index(ref_transcripts, index_folder, kmer_len=kmer_len)
    exec_this(cmds)


def do_quant_sailfish_single_end(index_dir, libtype, single_end_seq_read="", quant_dir="."):
    cmds = build_cmd_sailfish_quant(index_dir=index_dir, libtype=libtype, unmated_seq=single_end_seq_read, quant_dir=quant_dir)
    tmpp(" ".join(cmds))
    exec_this(cmds)


def do_quant_sailfish_pair_end(index_dir, libtype, mate1="", mate2="", quant_dir="."):
    cmds = build_cmd_sailfish_quant(index_dir, libtype, mate1=mate1, mate2=mate2, quant_dir=quant_dir)
    tmpp(" ".join(cmds))
    exec_this(cmds)


class PipeLine():
    def __init__(self, list_console_cmd):
        self._customized_parameter = functools.partial(parse_parameters,
                                                       short_option_definition="c:g:a:r:1:2:o:k:h",
                                                       long_option_definition=["libtype", "help"])

        if list_console_cmd:

            opts, args = self._customized_parameter(cmd_args=list_console_cmd)
            self._assign_args(opts)

    def _assign_args(self, opts):
        setting_map_of_opts = dict(opts)
        current_dir = os.path.abspath(os.curdir)

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

        elif "-1" in setting_map_of_opts and "-2" in setting_map_of_opts and setting_map_of_opts["-1"] and setting_map_of_opts["-2"]:
            # yes , this is paired end operation
            self._is_reads_paired_end = True
            self._paired_end_reads_1_file_path = setting_map_of_opts["-1"]
            self._paired_end_reads_2_file_path = setting_map_of_opts["-2"]
        else:   # neither single end nor paired end
            print("illegal specification of sequencing reads file, exiting .....")
            sys.exit(-1)

        self._kmer_len = 19
        if "-k" in setting_map_of_opts and setting_map_of_opts["-k"]:
            try:
                self._kmer_len = int(setting_map_of_opts["-k"])
            except ValueError as e:
                print("you should input a integer to 'k', exiting ......")
                sys.exit(-1)

        self._index_libtype = "IU"
        if "--libtype" in setting_map_of_opts and setting_map_of_opts["--libtype"]:
            self._index_libtype = setting_map_of_opts["--libtype"]

        self._is_circ_isoform_provided =False
        if "-c" in setting_map_of_opts and setting_map_of_opts["-c"]:
            self._is_circ_isoform_provided = True
            self._ciri_output = setting_map_of_opts["-c"]
        else:
            # no CIRI means you have to do it with only linear transcription
            pass

    def process_the_pipe_line(self):
        # this the main part of perform the pipeline.
        if self._is_circ_isoform_provided and self._ciri_output:
            tmpp("do circular run")
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
        circular_only_transcript_fasta_file_path = os.path.join(self._output_folder, "circular_only_transcript.fa")
        whole_transcript_file_path = os.path.join(self._output_folder, "transcriptome_with_circular.fa")

        sailfish_index_folder_circular = os.path.join(self._output_folder, "index_circular")
        sailfish_quant_folder_circular = os.path.join(self._output_folder, "quant_circular")

        if not os.path.exists(linear_transcript_fasta_file_path):
            tmpp("extract linear transcript ..")
            do_extract_classic_linear_transcript(self._annotation_file, self._genomic_seq, linear_transcript_fasta_file_path)
            do_convert_in_site(linear_transcript_fasta_file_path)

        # prepare circular transcriptome fa file.

        tmpp("make a gtf for circular RNA")
        do_make_gtf_for_circ_exon(self._annotation_file, self._ciri_output, circular_only_annotation_exon_file_path)

        tmpp("extract circular RNA transcriptome")
        do_extract_circular_transcript(circular_only_annotation_exon_file_path, self._genomic_seq, circular_only_transcript_fasta_file_path)
        do_convert_in_site(circular_only_transcript_fasta_file_path)
        do_add_addapt(circular_only_transcript_fasta_file_path)

        do_combine_circular_fa(linear_transcript_fasta_file_path, circular_only_transcript_fasta_file_path, whole_transcript_file_path)

        tmpp("make index for sailfish")
        do_make_index_sailfish(ref_transcripts=whole_transcript_file_path, index_folder=sailfish_index_folder_circular, kmer_len=self._kmer_len)

        if self._is_reads_paired_end:
            tmpp("do pair end reads quantification")
            do_quant_sailfish_pair_end(sailfish_index_folder_circular,
                                       libtype=self._index_libtype,
                                       mate1=self._paired_end_reads_1_file_path,
                                       mate2=self._paired_end_reads_2_file_path,
                                       quant_dir=sailfish_quant_folder_circular)

        else:
            tmpp("do single end reads quantification")
            do_quant_sailfish_single_end(sailfish_index_folder_circular,
                                         libtype=self._index_libtype,
                                         single_end_seq_read=self._single_end_read_file_path,
                                         quant_dir=sailfish_quant_folder_circular)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print(__doc__)
    else:
        work_on_it = PipeLine(sys.argv[1:])
        work_on_it.process_the_pipe_line()