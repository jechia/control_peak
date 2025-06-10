#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genomic Peak Analysis Tool

Analyzes overlap between genomic peaks and transcript annotations,
generating control regions for comparative analysis.

peak example
chr start stop dataset_label 1000 strand log2(eCLIP fold-enrichment over size-matched input) -log10(eCLIP vs size-matched input p-value) -1 -1 chr start end name 0 strand
chr7    151087809   151087876   EWSR1_K562_rep01    1000    +   4.49672532155361    16.2109119555925    -1  -1  chr7    151086742   151123635   AGAP3-211   0   +

annotation example
chr1    67092164    67093004    C1orf141-203_UTR3_67092164  0
@author: huyue
"""

from argparse import ArgumentParser
import time
import multiprocessing
from itertools import zip_longest, chain
import intervals as I
import random
from operator import itemgetter

# Global constants
REGIONS = ("UTR3", "UTR5", "CDS", "exon")
MIN_REGION_LENGTH = 50
DEFAULT_EXTENSION_LENGTH = 100
SLIDING_WINDOW_SIZE = 200
SLIDING_WINDOW_STEP = 50

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def create_interval(interval_str,is_list=True):
    """
    Create interval object from string or list of strings.
    
    Args:
        interval_str: String like "100-200" or list of such strings
        is_list: Whether input is a list
        
    Returns:
        intervals.Interval object
    """
    interval=I.empty()

    if is_list:
        for coord_str in interval_str:
            start, end = coord_str.split("-")
            interval = interval|I.closed(int(start),int(end))
    else:
        start, end = interval_str.split("-")
        interval=I.closed(int(start),int(end))
    
    return interval

def findoverlap(peak_range,annotation_info):
    """
    Find overlaps between peak and annotation regions.
    
    Args:
        peak_range: Peak coordinate string "start-end"
        annotation_info: Dict of region -> coordinate list
        
    Returns:
        Tuple of (coverage_list, coverage_ratios, region_list)
    """
    peak_interval = create_interval(peak_range,is_list=False)
    coverage_list = []
    coverage_ratios = []
    region_list = []

    for region in annotation_info:
        region_interval = create_interval(annotation_info[region])
        intersection = peak_interval & region_interval
        if intersection != I.empty():
            for interval_part in list(intersection):
                overlap_start, overlap_end = list((interval_part.lower,interval_part.upper))
                coverage_ratio = (overlap_end - overlap_start) /(peak_interval.upper-peak_interval.lower)
                if region == "intron":
                    coverage_list.append(str(peak_interval.lower)+"-"+str(peak_interval.upper))
                else:
                    coverage_list.append(str(overlap_start)+"-"+str(overlap_end))
                coverage_ratios.append(coverage_ratio)
                region_list.append(region)
    return coverage_list,coverage_ratios,region_list

def build_peak_dict(peak_info):
    """
    Build dictionary of peaks organized by chromosome.
    
    Args:
        peak_info: List of peak file lines
        
    Returns:
        Dict: {chromosome: {peak_range: {transcript: [names], FC: fold_change}}}
    """
    peak_dict={}

    for line in peak_info:
        fields = line.strip('\n').split('\t')
        
        chrom = fields[0]
        transcript_name = fields[11]
        strand = fields[5]
        fold_change = fields[6]

        if strand == "+":
            start=int(fields[1])
            peak_range=str(start)+"-"+fields[2]
        elif strand == "-":
            end=int(fields[2])
            peak_range=fields[1]+"-"+str(end)
        else:
            continue

        # Initialize nested dictionaries
        if chrom not in peak_dict:
            peak_dict[chrom] = {}
        if peak_range not in peak_dict[chrom]:
            peak_dict[chrom][peak_range]={
                "transcript": [],
                "FC": fold_change
            }

        peak_dict[chrom][peak_range]["transcript"].append(transcript_name)
    return peak_dict
                                
def build_annotation_dict(anno):
    """
    Build annotation dictionary from annotation file lines.
    
    Args:
        anno_lines: List of annotation file lines
        
    Returns:
        Dict: {chromosome: {transcript: {region: [coordinates], strand: strand, gene: gene_name}}}
    """
    anno_dict = {}
    for line in anno:
        fields = line.strip('\n').split('\t')
        
        chrom=fields[0]
        coord_range = fields[1]+"-"+fields[2]
        name_parts = fields[3].split('_')
        
        transcript_name=name_parts[0]
        region=name_parts[1]
        strand=fields[5]

        if chrom not in anno_dict:
            anno_dict[chrom]={}
        if transcript_name not in anno_dict[chrom]:
            anno_dict[chrom][transcript_name] = {}
        if region not in anno_dict[chrom][transcript_name]:
            anno_dict[chrom][transcript_name][region]=[]
        if "strand" not in anno_dict[chrom][transcript_name]:
            anno_dict[chrom][transcript_name]["strand"]=strand
        anno_dict[chrom][transcript_name][region].append(coord_range)
        anno_dict[chrom][transcript_name]["gene"]=fields[6]
    return anno_dict

def build_gene_dict(gene):
    """
    Build gene dictionary from gene file lines.
    
    Args:
        gene_lines: List of gene file lines
        
    Returns:
        Dict: {chromosome: {gene_name: coordinate_range}}
    """
    gene_dict={}

    for line in gene:
        fields = line.strip('\n').split('\t')
        
        chrom = fields[0]
        coord_range = fields[1]+"-"+fields[2]
        gene_name=fields[3]

        if chrom not in gene_dict:
            gene_dict[chrom]={}
        if gene_name not in gene_dict[chrom]:
            gene_dict[chrom][gene_name]={}
        gene_dict[chrom][gene_name]=coord_range
    return gene_dict

class GeneAnnotation:
    """Class for handling gene annotations and generating control regions."""
    
    def __init__(self,annotation_info):
        exon_list=[]
        if annotation_info.get("exon") is not None:
            exon_list=annotation_info["exon"]
        self.exons=exon_list
        self.strand=annotation_info['strand']
        self.annotation_dict = {}
        self.introns = []
        self.splice_sites = []
        self.start = 0
        self.end = 0

    def _get_splice_sites(self):
        """Calculate splice sites from exons."""
        splice_sites = []
        for exon_range in self.exons:
            start, end = exon_range.split("-")
            splice_sites.append(int(start))
            splice_sites.append(int(end))
        
        splice_sites.sort()
        self.start = splice_sites[0]
        self.end = splice_sites[-1]
        self.splice_sites=splice_sites
    
    def _get_introns(self):
        """Calculate intron coordinates from exons."""
        self._get_splice_sites()
        introns=[]
        
        # Create introns from internal splice sites
        middle_sites = self.splice_sites[1:-1]
        for start,end in grouper(middle_sites,2):
            introns.append(str(start)+"-"+str(end))
        self.introns=introns

    def add_introns(self,annotation_info):
        """Add intron information to the annotation."""
        self._get_introns()

        # Build comprehensive annotation dictionary
        for region in REGIONS:
            if annotation_info.get(region) is not None:
                self.annotation_dict[region] = annotation_info[region][:]
        
        self.annotation_dict["intron"] = self.introns[:]

    def peak_extension(self,peak,region_list,region):
        """Extend peak if necessary"""
        peak_coords = list((int(peak.lower),int(peak.upper)))
        peak_length = peak_coords[1]-peak_coords[0]
        
        for reg in list(region_list):
            if peak & reg != I.empty():
                intersection = peak & reg
                region_coords = list((reg.lower,reg.upper))
                intersect_coords = list((intersection.lower,intersection.upper))
                break

        region_length = region_coords[1]-region_coords[0]
        inter_length=intersect_coords[1]-intersect_coords[0]
        
        if region_length > MIN_REGION_LENGTH:
            if peak_length == inter_length:
                final_length = DEFAULT_EXTENSION_LENGTH if region_length > DEFAULT_EXTENSION_LENGTH else region_length
                extend_length = final_length - peak_length 

                peak_start = intersect_coords[0] - int(extend_length/2)
                new_start =  peak_start if peak_start > region_coords[0] else region_coords[0]
                peak_end = new_start + DEFAULT_EXTENSION_LENGTH 
                new_end = peak_end if peak_end < region_coords[1] else region_coords[1]

                if new_end - new_start < MIN_REGION_LENGTH:
                    peak_start = new_end-final_length
                    new_start = peak_start if peak_start > region_coords[0] else region_coords[0]
                new_peak = str(new_start)+"-"+str(new_end)
                new_peak = create_interval(new_peak,is_list=False)
                return new_peak
            else:
                if region == "intron":
                    extend_length = DEFAULT_EXTENSION_LENGTH - peak_length
                    new_start = intersect_coords[0] - int(extend_length/2)
                    new_end = new_start + DEFAULT_EXTENSION_LENGTH
                    new_peak = str(new_start)+"-"+str(new_end)
                    new_peak = create_interval(new_peak,is_list=False)
                    return new_peak
                else:
                    return None
    def get_control(self,peak,region,control_list):
        """
        Generate control region
        """

        region_list=create_interval(self.annotation_dict[region])
        peak_coords=list((int(peak.lower),int(peak.upper)))
        peak_length=peak_coords[1]-peak_coords[0]
        
        if peak_length < MIN_REGION_LENGTH:
            return None,None,None,None,None
        
        if control_list == I.empty() or peak_length <= 0:
            return None, None, None, None, None
        
        if region in ["UTR3", "UTR5", "exon"]:
            return self._get_simple_control(peak_coords, peak_length, control_list, region)

        elif region == "CDS":
            return self._get_cds_control(peak_coords, peak_length, control_list, region_list)
        
        elif region == "intron":
            return self._get_intron_control(peak_coords, peak_length, control_list, region_list)

        return None, None, None, None, None
    
    def _get_simple_control(self, peak_coords, peak_length, control_list, region):
        """Generate control for simple regions (UTR3, UTR5, exon)."""
        start_list = []
        for control_region in list(control_list):
            control_region = list((int(control_region.lower), int(control_region.upper)))
            region_length = control_region[1] - control_region[0]
            if region_length > peak_length:
                start = range(control_region[0], control_region[1] - peak_length)
                start_list = chain(start_list, start)
        start_list = list(start_list)
        if len(start_list) > 0:
            control_start = random.choice(start_list)
            return peak_coords[0], peak_coords[1], control_start, control_start + peak_length, region
        
        return None, None, None, None, None
    
    def _get_cds_control(self, peak_coords, peak_length, control_list, region_list):
        """Generate control for CDS regions with splice site awareness."""
        peak_interval = create_interval(f"{peak_coords[0]}-{peak_coords[1]}", is_list=False)

        for exon in list(region_list):
            if peak_interval & exon != I.empty():
                exon_coords = list((exon.lower, exon.upper))
        
        # Calculate distances to splice sites
        dist_to_start = peak_coords[0] - exon_coords[0]
        dist_to_end = exon_coords[1] - peak_coords[1]

        # Choose splice sites based on proximity
        if dist_to_start < dist_to_end:
            splice_sites = self.splice_sites[1::2] # 3' splice sites
            distance = dist_to_start
        else:
            splice_sites = self.splice_sites[::2] # 5' splice sites
            distance = dist_to_end

        # Generate control intervals at corresponding positions
        control_intervals = I.empty()
        for splice_site in splice_sites: 
            control_pos = splice_site + distance
            control_interval = I.closed(control_pos, control_pos + peak_length)

            if control_interval & control_list == control_interval: 
                control_intervals = control_intervals | control_interval
        if control_intervals != I.empty():
            chosen_control = random.choice(list(control_intervals))
            control_coords = list((int(chosen_control.lower), int(chosen_control.upper)))
            return peak_coords[0], peak_coords[1], control_coords[0], control_coords[1], "CDS"
        
        return None, None, None, None, None
    
    def _get_intron_control(self, peak_coords, peak_length, control_list, region_list):
        """Generate control for intron regions"""

        control_intervals = I.empty()
        intersect_length = 0
        is_exon_overlap = False
        is_5_prime_ss = False

        # Check overlap with intron regions
        peak_interval = create_interval(f"{peak_coords[0]}-{peak_coords[1]}", is_list=False)
        for intron in list(region_list):
            if peak_interval & intron != I.empty():
                intersection = peak_interval & intron
                region_coords = list((intron.lower, intron.upper))
                intersect_coords = list((intersection.lower, intersection.upper))
                intersect_length = intersect_coords[1] - intersect_coords[0]

        # If no intron overlap, check exon overlap
        if intersect_length == 0:
            exon_list = create_interval(self.annotation_dict["exon"])
            for exon in list(exon_list):
                if peak_interval & exon != I.empty():
                    intersection = peak_interval & exon
                    region_coords = [exon.lower, exon.upper]
                    intersect_coords = [intersection.lower, intersection.upper]
                    intersect_length = intersect_coords[1] - intersect_coords[0]
                    is_exon_overlap = True
        
        if intersect_length == 0:
            return None, None, None, None, None
        
        # Calculate distances to region boundaries
        dist_to_start = peak_coords[0] - region_coords[0]
        dist_to_end = region_coords[1] - peak_coords[1]

        # Determine region type and splice sites
        if peak_length == intersect_length:
            if dist_to_start < dist_to_end:
                if is_exon_overlap:
                    splice_sites = self.splice_sites[::2]
                    final_region = "3pexon" if self.strand == "+" else "5pexon"
                    distance = dist_to_start
                else:
                    splice_sites = self.splice_sites[1::2]
                    final_region = "5pintron" if self.strand == "+" else "3pintron"
                    distance = dist_to_start
            else:
                if is_exon_overlap:
                    splice_sites = self.splice_sites[1::2]
                    final_region = "5pexon" if self.strand == "+" else "3pexon"
                    distance = dist_to_end
                    is_5_prime_ss = True
                else:
                    splice_sites = self.splice_sites[2::2]
                    final_region = "3pintron" if self.strand == "+" else "5pintron"
                    distance = dist_to_end
                    is_5_prime_ss = True
        else:
            if dist_to_start < 0:
                splice_sites = self.splice_sites[1::2]
                final_region = "5pss" if self.strand == "+" else "3pss"
                distance = dist_to_start
            if dist_to_end < 0:
                splice_sites = self.splice_sites[2::2]
                final_region = "3pss" if self.strand == "+" else "5pss"
                distance = dist_to_end
                is_5_prime_ss = True
        
        # Generate control intervals based on splice sites
        for splice_site in splice_sites:
            if is_5_prime_ss:
                control_pos = splice_site - distance
                control_interval = I.closed(control_pos - peak_length, control_pos)
            else:
                control_pos = splice_site + distance
                control_interval = I.closed(control_pos, control_pos + peak_length)
        
            if control_interval & control_list == control_interval:
                control_intervals = control_intervals | control_interval

        if control_intervals != I.empty():
            chosen_control = random.choice(list(control_intervals))
            control_coords = list((int(chosen_control.lower), int(chosen_control.upper)))
            return peak_coords[0], peak_coords[1], control_coords[0], control_coords[1], final_region
        
        return None, None, None, None, None


def sliding_window(sequence, window_size, min_len):
    """Generate sliding windows over a sequence with minimum length check."""
    result = []
    seq_len = len(sequence)

    for i in range(0, seq_len, window_size):
        if i + window_size > seq_len:
            if seq_len - i > min_len:
                window_str = str(sequence[i])+"-"+str(sequence[seq_len-1])
        else:
            window_str = str(sequence[i])+"-"+str(sequence[i+window_size])
        result.append(window_str)

    return result

def parsePeaks(peak_chr_dict,anno_chr_dict,gene_chr_dict,chromosome):
    """
    Process peaks for a single chromosome
    """

    results = []
    processed_peaks_dict = {}
    coverage_dict = {}
    assemble_dict = {}
    peak_region_dict = {}

    total_peak_dict = {
        "list": {"plus": [], "minus": []},
        "interval": {"plus": I.empty(), "minus": I.empty()}
    }

    total_control_dict = {
        "list": {"plus": [], "minus": []},
        "interval": {"plus": I.empty(), "minus": I.empty()}
    }

    # First pass: build assemble and coverage dict
    for peak_range in peak_chr_dict:
        peak_info = peak_chr_dict[peak_range]

        for transcript_name in peak_info["transcript"]:
            if transcript_name not in anno_chr_dict:
                continue

            annotatoin = anno_chr_dict[transcript_name]
            gene_anno = GeneAnnotation(annotatoin)
            gene_anno.add_introns(annotatoin)

            strand_key = "plus" if annotatoin['strand'] == "+" else "minus"
            total_peak_dict["list"][strand_key].append(peak_range)

            overlaps, coverage_ratios, region_types = findoverlap(peak_range, gene_anno.annotation_dict)
            if len(coverage_ratios) > 0:
                for i, (overlap, coverage, region_type) in enumerate(zip(overlaps, coverage_ratios, region_types)):
                    # BUild coverage dict
                    if transcript_name not in coverage_dict:
                        coverage_dict[transcript_name] = {"all": 0}
                    if region_type not in coverage_dict[transcript_name]:
                        coverage_dict[transcript_name][region_type] = 0
                    peak_interval = create_interval(peak_range, is_list=False)
                    peak_length = peak_interval.upper - peak_interval.lower

                    coverage_dict[transcript_name]["all"] += peak_length
                    coverage_dict[transcript_name][region_type] += coverage * peak_length

                    # Build assemble dict
                    if transcript_name not in assemble_dict:
                        assemble_dict[transcript_name] = {"all": []}
                    if region_type not in assemble_dict[transcript_name]:
                        assemble_dict[transcript_name][region_type] = []
                    
                    assemble_dict[transcript_name][region_type].append(peak_range)
                    assemble_dict[transcript_name]["all"].append(peak_range)

                    # Track peaks with complete overlap or intron peaks
                    if coverage == 1 or region_type == "intron":
                        if peak_range not in peak_region_dict:
                            peak_region_dict[peak_range] = {}
                        if region_type not in peak_region_dict[peak_range]:
                            peak_region_dict[peak_range][region_type] = []
                        peak_region_dict[peak_range][region_type].append(transcript_name)
    
    # Build peak interval sets
    total_peak_dict["interval"]["plus"] = create_interval(total_peak_dict["list"]["plus"])
    total_peak_dict["interval"]["minus"] = create_interval(total_peak_dict["list"]["minus"])

    # Second pass: generate controls for qualified peaks
    for peak_range in peak_region_dict:
        # Determine peak region priority
        if "CDS" in peak_region_dict[peak_range]:
            peak_region = "CDS"
        elif "UTR3" in peak_region_dict[peak_range]:
            peak_region = "UTR3"
        elif "UTR5" in peak_region_dict[peak_range]:
            peak_region = "UTR5"
        elif "exon" in peak_region_dict[peak_range]:
            peak_region = "exon"
        else:
            peak_region = "intron"

        peak_interval = create_interval(peak_range, is_list=False)
        peak_fc=peak_chr_dict[peak_range]["FC"]

        for transcript_name in peak_region_dict[peak_range][peak_region]:
            if transcript_name not in anno_chr_dict:
                continue
            
            annotation = anno_chr_dict[transcript_name]
            strand_key = "plus" if annotation['strand'] == "+" else "minus"

            gene_anno = GeneAnnotation(annotation)
            gene_anno.add_introns(annotation)
            gene_name = annotation["gene"]

            # Determine available control regions\
            if gene_name in gene_chr_dict:
                gene_bounds = create_interval(gene_chr_dict[gene_name], is_list=False)
                bound_list = gene_bounds & total_peak_dict["interval"][strand_key]
            else:
                bound_list = total_peak_dict["interval"][strand_key]
            
            # Get region list for control generation
            if peak_region == "intron":
                region_interval = create_interval(f"{gene_anno.start}-{gene_anno.end}", is_list=False)
            else:
                region_interval = create_interval(gene_anno.annotation_dict[peak_region])
            
            control_list = region_interval - bound_list
            peak_coords = list((int(peak_interval.lower), int(peak_interval.upper)))
            peak_length = peak_coords[1] - peak_coords[0]

            # Process peak - either as whole or in sliding windows
            peak_processed = False
            windows = ([f"{peak_coords[0]}-{peak_coords[1]}"] if peak_length <= 400 
                      else sliding_window(range(peak_coords[0], peak_coords[1] + 1), 
                                       SLIDING_WINDOW_SIZE, SLIDING_WINDOW_STEP))
            
            for window in windows:
                window_interval = create_interval(window, is_list=False)
                if window_interval in processed_peaks_dict:
                    continue

                peak_start, peak_end, control_start, control_end, region_type = gene_anno.get_control(window_interval, peak_region, control_list)

                if not control_start:
                    continue

                control_interval = create_interval(f"{control_start}-{control_end}", is_list=False)
                control_overlap1 = control_interval & total_peak_dict["interval"][strand_key]
                control_overlap2 = control_interval & total_control_dict["interval"][strand_key]


                if control_overlap1 == I.empty() and control_overlap2 == I.empty() and window_interval not in processed_peaks_dict:
                    control_range = f"{control_start}-{control_end}"
                    total_control_dict["list"][strand_key].append(control_range)
                    total_control_dict["interval"][strand_key] = create_interval(total_control_dict["list"][strand_key])

                    processed_peaks_dict[window_interval] = control_start
                    peak_name = f"{gene_name}_{region_type}_{peak_fc}"
                    result_fields = [
                        chromosome,
                        str(peak_start),
                        str(peak_end),
                        peak_name,
                        annotation['strand'],
                        str(control_start),
                        str(control_end)
                    ]
                    results.append("\t".join(result_fields))
                    peak_processed = True
            
            if peak_processed:
                break

    return results
               
def createArgs():
    """
    Create the command line interface of the program.
    """
    parser = ArgumentParser(description="Analyze genomic peaks and generate controls")  
    parser.add_argument("-i","--input",dest="input",help="Input peaks file",required=True)
    parser.add_argument("-a","--anno",dest="anno",help="Transcript annotation file",required=True)
    parser.add_argument("-g","--gene",dest="gene",help="Gene annotation file",required=True)
    parser.add_argument("-p","--pool",dest="pool",type=int,help="Number of cores for multiprocessing",required=True)
    parser.add_argument("-o", "--output", dest="output", help="Output directory", default="./new_result/")
    parser.add_argument("-s", "--seed", dest="seed", help="Seed for random selection", type=int,default=1234)
    return parser
    return parser

if __name__ == '__main__':
    start_time=time.time()

    # Parse arguments
    args = createArgs().parse_args()

    random.seed(args.seed)

    print('Building annotation dict...')

    # Load annoation files
    with open(args.anno) as f:
        anno_info = f.readlines()
    anno_dict = build_annotation_dict(anno_info)
    print(f"Annotation dictionary built: {time.time() - start_time:.2f} s")

    with open(args.gene) as f:
        gene_info = f.readlines()
    gene_dict = build_gene_dict(gene_info)
    print(f"Gene dictionary built: {time.time() - start_time:.2f} s")

    with open(args.input) as f:
        peak_info = f.readlines()
    peak_dict = build_peak_dict(peak_info)
    print(f"Peak dictionary built: {time.time() - start_time:.2f} s")

    pool = multiprocessing.Pool(processes=args.pool)
    total_results = []

    # start processing peaks

    for chromosome in peak_dict:
        if chromosome in anno_dict:
            peak_chr_dict = peak_dict[chromosome]
            anno_chr_dict = anno_dict[chromosome]
            gene_chr_dict = gene_dict.get(chromosome, {})
            result = pool.apply_async(parsePeaks, args=(peak_chr_dict, anno_chr_dict, gene_chr_dict, chromosome))
            total_results.append(result)
    
    pool.close()
    pool.join()
    print(f"Peak processing completed: {time.time() - start_time:.2f} s")

    # Write output
    import os
    os.makedirs(args.output, exist_ok=True)
    
    base_name = os.path.splitext(os.path.basename(args.input))[0]
    output_file = os.path.join(args.output, f"{base_name}_control.bed")

    with open(output_file, "w") as fnOut:
        for result in total_results:
            event_list = result.get()
            for event in event_list:
                fnOut.write(event + '\n')

    print(f'Analysis completed successfully!')
    print(f'Output written to: {output_file}')
    print(f'Total time: {time.time() - start_time:.2f}s')
    

