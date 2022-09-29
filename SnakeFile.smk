import os.path
import pysam
import numpy as np
import pandas as pd


configfile: "config.yaml"
version = "0.1.Alpha"
print("MERD: Meiotic Recombination Detection Pipeline with PacBio HiFi sequencing reads")
ava_family = {}
unavailable_family = []
available_family = {}
targets = []
target = "{family}.{sample}.{hap}.{tech}"
for family_name, info in config["family"].items():
    print(f"------------------------------")
    print(f"[Info] Checking {family_name}")

    if ("merged_vcf" not in info) or (not os.path.exists(info["merged_vcf"])):
        unavailable_family.append(family_name)
        print(f"[Err] vcf file {info['merged_vcf']} not available! ")
        continue
    else:
        vcf_samples = [s for s in pysam.VariantFile(info["merged_vcf"]).header.samples]
        print(f"[Info] Samples in vcf file: {', '.join(vcf_samples)}")

    available_bam_samples = {}
    # print(info)
    for sample, tech_info in info["bams"].items():
        sample_bam = []
        for tech, bam in tech_info.items():
            if not os.path.exists(bam):
                print(f"[Warning] bam file {bam} not available! ")
            else:
                sample_bam.append(tech)
        if len(sample_bam) > 0:
            available_bam_samples[sample] = sample_bam
    available_hap = []
    for hap in ["maternal", "paternal"]:
        if hap not in info["names"]:
            continue
            print(f"[War] no {hap} state in config.yaml file!")
        else:
            sample_name = info['names'][hap] if isinstance(info['names'][hap],str) else list(info['names'][hap].keys())[
                0]

            if sample_name not in vcf_samples:
                print(f"[War] {hap} sample {sample_name} not available  in vcf file {info['merged_vcf']}")
                continue
            elif sample_name not in available_bam_samples:
                print(f"[War] no available bam filr for {hap} sample {sample_name} ")
                continue
            else:
                available_hap.append(hap)
                print(f"[Info] {hap}: {sample_name}")
    for proband in info["names"]["probands"]:
        proband_sample_name = proband if isinstance(proband,str) else list(proband.keys())[0]
        if proband_sample_name not in vcf_samples:
            print(f"[War] proband sample {proband_sample_name} not available  in vcf file {info['merged_vcf']}")
            continue
        else:
            print(f"[Info] proband: {proband_sample_name}")
            for hap in available_hap:
                hap_name = info['names'][hap] if isinstance(info['names'][hap],str) else list(info['names'][hap])[0]
                techs = available_bam_samples[hap_name]
                print(f"[Info]     \t{hap} tech: {', '.join(techs)}")
                for tech in techs:
                    if family_name not in available_family:
                        available_family[family_name] = {}
                    if proband_sample_name not in available_family[family_name]:
                        available_family[family_name][proband_sample_name] = {}
                    if hap not in available_family[family_name][proband_sample_name]:
                        available_family[family_name][proband_sample_name][hap] = []
                    available_family[family_name][proband_sample_name][hap].append(tech)
                    if tech == "ONT": continue
                    targets.append(f"{family_name}/{family_name}.{proband_sample_name}.{hap}.{tech}")
# print("\n".join(targets))
dir_work = config["dir_work"] if "dir_work" in config else "./"
chromsomes = config["chromsomes"].split(",")
path_ref = config["reference"]
whatshap = config["software"]["whatshap"] if "whatshap" in config["software"] else "whatshap"
tabix = config["software"]["tabix"] if "tabix" in config["software"] else "tabix"
wildcard_constraints:
    chrom="|".join(chromsomes),
    hap="paternal|maternal",
    tech="HiFi|ONT|ONTUL"

rule all:
    input:
        expand(dir_work + "{target}.final.merd.tsv",target=targets)


# dir_work+ "{family}/{chrom}/{family}.{proband}.{hap}.{tech}.{chrom}.phased.byreads.csv",


# expand(dir_work + "{family}/{family}.ped",family="ChineseQuartet")


def get_merged_vcf(wildcards):
    vcf = config["family"][wildcards.family]["merged_vcf"]
    return {"vcf": vcf, "vcf_idx": vcf + ".tbi"}


rule filter_vcf:
    input:
        unpack(get_merged_vcf)
    output:
        vcf=dir_work + "{family}/{family}.pass.vcf.gz",
    run:
        vcf_input = pysam.VariantFile(f"{input.vcf}")
        samples = vcf_input.header.samples
        family_name = wildcards.family
        samples_dict = config["family"][wildcards.family]["names"]
        # print(samples_dict)
        hap_info = {}
        info = config["family"][wildcards.family]
        for hap in ["maternal", "paternal"]:
            if hap not in config["family"][wildcards.family]["names"]:
                # hap_info[hap] = "."
                pass
            else:
                sample_name = info['names'][hap] if isinstance(info['names'][hap],str) else \
                    list(info['names'][hap].keys())[0]
                if sample_name not in samples: continue
                hap_info[hap] = sample_name
        proband_names = [list(proband.keys())[0] if (not isinstance(proband,str)) else proband for proband in
                         samples_dict["probands"]]
        par_sample = hap_info["paternal"] if "paternal" in hap_info else None
        mar_sample = hap_info["maternal"] if "maternal" in hap_info else None

        if par_sample == None and mar_sample == None:
            print("[Err] No available parental samples!")
            exit()
        elif par_sample == None or mar_sample == None:
            parent_sample = mar_sample if par_sample == None else par_sample
            parent_num = 1
        else:
            parent_num = 2
        chroms = chromsomes
        # chrom_lens = {contig: info.length for contig, info in vcf_input.header.contigs.items() if contig in chroms}
        vcf_output_pass = pysam.VariantFile(str(output.vcf),"w",header=vcf_input.header)
        af_thresh = 0.2
        af_thresh_low = 0.05
        num_somatic = 0
        num_denovel = 0
        num_pass = 0
        num_pass2 = 0
        depth_dict = {sample: [] for sample in samples}
        for rec in vcf_input.fetch(chroms[0]):
            for sample, info in rec.samples.items():
                if sample in depth_dict and info["DP"] != None:
                    depth_dict[sample].append(info["DP"])
        depth_value = {sample: np.median(dps) for sample, dps in depth_dict.items()}
        num = 0
        for rec in vcf_input.fetch():
            num += 1
            if rec.contig not in chroms: continue  # 过滤非染色体
            ref_str = set(rec.ref)
            if len([j for i in rec.alts for j in i]) > 1: continue  # 过滤多个基因型的，只保留两种基因型的
            gt_dict = {}
            dp_dict = {}
            filter = False
            for sample, info in rec.samples.items():
                if None == info["DP"] or None in info["GT"]:
                    filter = True
                    continue
                else:
                    if info["DP"] > 3 * depth_value[sample] or info["DP"] < depth_value[sample] * 0.2:
                        filter = True
                        continue
                gt_dict[sample] = info["GT"]
                dp_dict[sample] = info["DP"]
            if filter:
                continue  # 过滤没有GT和DP不符合要求的记录

            # for child in proband_names:
            filter_genetic = False
            if parent_num == 2:
                for child in proband_names:
                    if not ((gt_dict[child][0] in gt_dict[par_sample] and gt_dict[child][1] in gt_dict[mar_sample]) or \
                            (gt_dict[child][1] in gt_dict[par_sample] and gt_dict[child][0] in gt_dict[mar_sample
                            ])):  # de nove mutation
                        filter_genetic = True
                        continue
            else:
                for child in proband_names:
                    if (gt_dict[child][0] not in gt_dict[parent_sample]) and (
                            gt_dict[child][1] not in gt_dict[parent_sample]):
                        filter_genetic = True
                        continue

            if filter_genetic: continue
            if wildcards.family == "ChineseQuartet":
                if not (((gt_dict["LCL5"][0] == gt_dict["LCL6"][0]) and (gt_dict["LCL5"][1] == gt_dict["LCL6"][1])) or \
                        ((gt_dict["LCL5"][0] == gt_dict["LCL6"][1]) and (
                                gt_dict["LCL5"][1] == gt_dict["LCL6"][0]))):  # twins 基因型不一致
                    continue
            vcf_output_pass.write(rec)
        vcf_output_pass.close()

rule get_ped:
    output:
        ped=dir_work + "{family}/{family}.ped",
    run:
        family_name = wildcards.family
        samples_dict = config["family"][wildcards.family]["names"]
        # print(samples_dict)
        output_ped = open(f"{output.ped}","w")
        hap_info = {}
        info = config["family"][wildcards.family]
        for hap in ["maternal", "paternal"]:
            if hap not in config["family"][wildcards.family]["names"]:
                hap_info[hap] = "."
            else:
                sample_name = info['names'][hap] if isinstance(info['names'][hap],str) else \
                    list(info['names'][hap].keys())[0]
                hap_info[hap] = sample_name
        for proband in samples_dict["probands"]:
            if (not isinstance(proband,str)):
                sample_name = list(proband.keys())[0]
                if "sex" in proband[sample_name]:
                    sex_id = "2" if proband[sample_name]["sex"] in ["Female", "FeMale", "F"] else "1"
                else:
                    sex_id = "1"
            else:
                sample_name = proband
                sex_id = "1"
            output_ped.write("\t".join([wildcards.family, sample_name, hap_info["paternal"], hap_info["maternal"],
                                        sex_id, "1"]) + "\n")
        output_ped.close()


# if proband[0]["sex"] in ["Female","F"] :
#     print("gkkgkgk")
# else:
#     print("lfllf")


def phasing_vcf_input(wildcards):
    vcf = dir_work + f"{wildcards.family}/{wildcards.family}.pass.vcf.gz"
    sample = config["family"][wildcards.family]["names"][wildcards.hap]
    sample = sample if isinstance(sample,str) else list(sample.keys())[0]
    bam = config["family"][wildcards.family]["bams"][sample][wildcards.tech]
    return {"vcf": vcf, "bam": bam, "bai": f"{bam}.bai", "vcf_idx": vcf + ".tbi"}


def get_bam_from_proband(wildcards):
    return config["family"][wildcards.family]["bams"][wildcards.proband][wildcards.tech]


rule phaing_vcf_ped:
    input:
        ped=dir_work + "{family}/{family}.ped",
        vcf=dir_work + "{family}/{family}.pass.vcf.gz",
        bam=get_bam_from_proband,
    output:
        dir_work + "{family}/{family}.{proband}.{tech}.ped.phased.vcf.gz",
    run:
        shell("{whatshap} phase  --indels -r {path_ref} --ped {input.ped}  -o {output} {input.vcf} {input.bam}")


rule phasing_vcf:
    input:
        unpack(phasing_vcf_input)
    output:
        vcf=temp(dir_work + "{family}/{chrom}/{family}.{proband}.{hap}.{tech}.{chrom}.phased.vcf.gz",)
    run:
        hap = wildcards.hap
        sample = config["family"][wildcards.family]["names"][wildcards.hap]
        sample = sample if isinstance(sample,str) else list(sample.keys())[0]
        shell("{whatshap} phase --chromosome {wildcards.chrom} --merge-reads --algorithm hapchat "
              "--ignore-read-groups --indels -r {path_ref} --sample {sample} -o {output} {input.vcf} {bam}")

rule phaing_vcf_ped_filter:
    input:
        vcf=dir_work + "{family}/{family}.{proband}.{tech}.ped.phased.vcf.gz",
    output:
        vcf=dir_work + "{family}/{family}.{proband}.{tech}.{hap}.ped.phased.vcf.gz",
    run:
        sample = config["family"][wildcards.family]["names"][wildcards.hap]
        sample = sample if isinstance(sample,str) else list(sample.keys())[0]
        print(sample)

        vcf = pysam.VariantFile(f"{input.vcf}","r")
        vcf_new = pysam.VariantFile(f"{output.vcf}","w",header=vcf.header)
        for rec in vcf.fetch():

            if (not rec.samples[sample].phased) or len(rec.alts) > 1: continue
            if abs(len(rec.ref) - len(rec.alts[0]) > 0): continue
            gt = rec.samples[sample]["GT"]
            if gt[0] == gt[1]: continue
            vcf_new.write(rec)
        vcf.close()
        vcf_new.close()


def get_bam(wildcards):
    sample = config["family"][wildcards.family]["names"][wildcards.hap]
    sample = sample if isinstance(sample,str) else list(sample.keys())[0]
    return config["family"][wildcards.family]["bams"][sample][wildcards.tech]


def get_dis(x, y):
    size = 0
    for i, j in zip(x,y):
        size += (abs(i - j))
    return size / len(x)


rule extract_shift_variants_from_HiFi_reads:
    input:
        vcf=dir_work + "{family}/{family}.{proband}.{tech}.{hap}.ped.phased.vcf.gz",
        vcf_idx=dir_work + "{family}/{family}.{proband}.{tech}.{hap}.ped.phased.vcf.gz.tbi",
        ref=config["reference"],
        bam=get_bam
    output:
        bam=dir_work + "{family}/details/{family}.{proband}.{hap}.{tech}.phased.byreads.{chrom}.bam",
        read_list=dir_work + "{family}/details/{family}.{proband}.{hap}.{tech}.phased.byreads.{chrom}.list",
    run:
        mapping_qual_thresh = config["mapping_qual_thresh"] if "mapping_qual_thresh" in config else 30
        sample = config["family"][wildcards.family]["names"][wildcards.hap]
        sample = sample if isinstance(sample,str) else list(sample.keys())[0]
        # sample = wildcards.proband
        vcf = pysam.VariantFile(f"{input.vcf}","r")
        bam = pysam.Samfile(f"{input.bam}","r")
        new_bam = pysam.Samfile(f"{output.bam}","w",header=bam.header)
        new_list = open(f"{output.read_list}","w")
        new_list.write("read.query_name\ttotal_mut_num\thap1_dis\thap2_dis\t"
                       "chrom\tread.reference_start\tread.reference_end\t"
                       f"read.mapping_quality\tread.is_supplementary\tread.is_secondary\tsour_vec\thap1_str\thap2_str\tread_str\tshift_pos\n")
        ref = pysam.FastaFile(f"{input.ref}")
        chroms = [f"chr{i}" for i in range(1,23)] + ["chrX"]
        chrom = wildcards.chrom
        chrom_lens = {contig: info.length for contig, info in vcf.header.contigs.items() if contig in chroms}
        # for chrom, chrom_len in chrom_lens.items():
        #     print(chrom)
        chrom_len = chrom_lens[chrom]
        mut_pots = np.zeros(chrom_len + 1)
        hap1_pots = np.zeros(chrom_len + 1)
        hap2_pots = np.zeros(chrom_len + 1)
        for rec in vcf.fetch(chrom):
            # if (not rec.samples[sample].phased) or len(rec.alts) > 1: continue
            gt = rec.samples[sample]["GT"]
            # if gt[0] == gt[1]: continue
            if gt[0] == 1:
                hap1_pots[rec.pos] = 1
            else:
                hap2_pots[rec.pos] = 1
            mut_pots[rec.pos] = 1
        # print(mut_pots.sum(),hap1_pots.sum(),hap2_pots.sum())
        for read in bam.fetch(chrom):
            # print(read.reference_start,read.reference_end)
            if read.mapping_quality < mapping_qual_thresh or read.is_secondary or read.is_duplicate: continue
            if mut_pots[read.reference_start:read.reference_end].sum() < 2: continue
            ref_muts_pos = [rec.pos for rec in vcf.fetch(chrom,read.reference_start,read.reference_end)]
            pairs = {item[1]: item[0] for item in read.get_aligned_pairs(matches_only=True)}
            total_mut_num = len(ref_muts_pos)
            vec_reads = []
            vec_hap1 = []
            vec_hap2 = []
            # print(ref_muts_pos)
            mut_pos_in = []
            for pos in ref_muts_pos:

                if pos not in pairs:
                    # varvec.append("N")
                    continue
                else:
                    mut_pos_in.append(pos)
                    if (hap1_pots[pos] == 0) and (hap2_pots[pos] == 0):
                        print("error position")
                        print(hap1_pots[pos - 2:pos + 2])
                        print(hap2_pots[pos - 2:pos + 2])
                        exit()
                    read_pos = pairs[pos]
                    read_base = read.query_sequence[read_pos - 1]
                    ref_base = ref.fetch(chrom,pos - 1,pos)
                    if read_base != ref_base:
                        vec_reads.append(1)
                    else:
                        vec_reads.append(0)
                    if hap1_pots[pos] == 1:
                        vec_hap1.append(1)
                        vec_hap2.append(0)
                    else:
                        vec_hap1.append(0)
                        vec_hap2.append(1)
            if len(mut_pos_in) < 2: continue
            hap1_dis = get_dis(vec_hap1,vec_reads)
            hap2_dis = get_dis(vec_hap2,vec_reads)
            sour_vec = ""
            for i, j in zip(vec_hap1,vec_reads):
                if i == j:
                    sour_vec += "1"
                else:
                    sour_vec += "2"
            vec_hap1_str = "".join([f"{i}" for i in vec_hap1])
            vec_hap2_str = "".join([f"{i}" for i in vec_hap2])
            vec_reads_str = "".join([f"{i}" for i in vec_reads])
            frez = sour_vec[0]
            last_i = sour_vec[0]
            shift_index = 0
            for i, index in zip(sour_vec,range(len(sour_vec))):
                if i == last_i:
                    continue
                else:
                    shift_index = index
                    frez += i
                    last_i = i
            if len(frez) != 2: continue
            # print(frez,sour_vec)
            if 0.0 < hap1_dis < 1:
                # print(read.query_name)
                new_bam.write(read)
                new_list.write(f"{read.query_name}\t{total_mut_num}\t{hap1_dis}\t{hap2_dis}\t"
                               f"{chrom}\t{read.reference_start}\t{read.reference_end}\t"
                               f"{read.mapping_quality}\t{read.is_supplementary}\t{read.is_secondary}\t{sour_vec}\t"
                               f"{vec_hap1_str}\t{vec_hap2_str}\t{vec_reads_str}\t{mut_pos_in[shift_index - 1]},{mut_pos_in[shift_index]}\n")
        new_bam.close()
        new_list.close()
        bam.close()
        vcf.close()
        ref.close()


def get_abnormal_filtering_input(wildcards):
    pat_sample = config["family"][wildcards.family]["names"]["paternal"]
    pat_sample = pat_sample if isinstance(pat_sample,str) else list(pat_sample.keys())[0]
    mat_sample = config["family"][wildcards.family]["names"]["maternal"]
    mat_sample = mat_sample if isinstance(mat_sample,str) else list(mat_sample.keys())[0]
    sample = config["family"][wildcards.family]["names"][wildcards.hap]
    sample = sample if isinstance(sample,str) else list(sample.keys())[0]
    bam = config["family"][wildcards.family]["bams"][sample][wildcards.tech]
    proband_bam = config["family"][wildcards.family]["bams"][wildcards.proband][wildcards.tech]
    mat_bam = config["family"][wildcards.family]["bams"][mat_sample][wildcards.tech]
    pat_bam = config["family"][wildcards.family]["bams"][pat_sample][wildcards.tech]
    read_list = dir_work + f"{wildcards.family}/details/{wildcards.family}.{wildcards.proband}.{wildcards.hap}.{wildcards.tech}.phased.byreads.{wildcards.chrom}.list",
    vcf = dir_work + f"{wildcards.family}/{wildcards.family}.{wildcards.proband}.{wildcards.tech}.{wildcards.hap}.ped.phased.vcf.gz",
    return {"pat_bam": pat_bam, "mat_bam": mat_bam, "bam": bam, "read_list": read_list, "vcf": vcf,
            "ref": config["reference"], "proband_bam": proband_bam}


def get_max2_gt(gt_dict, sample_depth):
    this_depth = np.sum([j for i, j in gt_dict.items()])
    if this_depth < sample_depth / 3:
        return [], this_depth

    sorted_gts = sorted(gt_dict.items(),key=lambda x: x[1],reverse=True)
    if len(sorted_gts) == 1:
        return [sorted_gts[0][0]] * 2, this_depth
    else:
        if sorted_gts[0][1] > this_depth * 0.75:
            return [sorted_gts[0][0]] * 2, this_depth
        elif sorted_gts[0][1] + sorted_gts[1][1] > this_depth * 0.75:
            return [sorted_gts[0][0], sorted_gts[1][0]], this_depth
        else:
            return [], this_depth


rule meiotic_recombination_detection:
    input:
        unpack(get_abnormal_filtering_input)
    output:
        read_list=dir_work + "{family}/details/{family}.{proband}.{hap}.{tech}.phased.byreads.pass.{chrom}.list",
    run:
        mapping_qual_thresh = config["mapping_qual_thresh"] if "mapping_qual_thresh" in config else 30
        chrom = wildcards.chrom
        sample = config["family"][wildcards.family]["names"][wildcards.hap]
        sample = sample if isinstance(sample,str) else list(sample.keys())[0]
        ref = pysam.FastaFile(f"{input.ref}")
        bam = pysam.Samfile(f"{input.bam}")
        pat_bam = pysam.Samfile(f"{input.pat_bam}")
        mat_bam = pysam.Samfile(f"{input.mat_bam}")
        proband_bam = pysam.Samfile(f"{input.proband_bam}")
        df_input = pd.read_csv(f"{input.read_list}",sep="\t")
        vcf = pysam.VariantFile(f"{input.vcf}")
        chroms = [f"chr{i}" for i in range(1,23)] + ["chrX"]

        chrom_lens = {contig: info.length for contig, info in vcf.header.contigs.items() if contig in chroms}
        chrom_len = chrom_lens[chrom]

        mut_pots = np.zeros(chrom_len + 1)
        hap1_pots = np.zeros(chrom_len + 1)
        hap2_pots = np.zeros(chrom_len + 1)
        for rec in vcf.fetch(chrom):
            # if (not rec.samples[sample].phased) or len(rec.alts) > 1: continue
            gt = rec.samples[sample]["GT"]
            # if gt[0] == gt[1]: continue
            if gt[0] == 1:
                hap1_pots[rec.pos] = 1
            else:
                hap2_pots[rec.pos] = 1
            mut_pots[rec.pos] = 1

        sample_depth_dict = {}
        for this_bam, bam_name in zip([bam, proband_bam, pat_bam, mat_bam],["this_hap", "proband", "pat", "mat"]):
            depths = []
            np.random.seed(100)
            for item in range(100):
                region_start = np.random.randint(0,chrom_len - 1000)
                region_end = region_start + 500
                this_dps = [i.nsegments for i in
                            bam.pileup(contig=chrom,start=region_start,end=region_end,
                                min_mapping_quality=mapping_qual_thresh,truncate=True)]
                if len(this_dps) == 0: continue
                this_region_depth = np.nanmean(this_dps)
                depths.append(this_region_depth)
            read_depth_mean = np.nanmean(depths)
            read_depth_std = np.nanstd(depths)
            sample_depth_dict[bam_name] = [read_depth_mean, read_depth_std]
        output_dict = {}
        for read_id, info in df_input.iterrows():
            if info["read.is_secondary"]: continue
            start, end = info["shift_pos"].split(",")
            id_index = "_".join([start, end])
            if id_index not in output_dict:
                output_dict[id_index] = []
            output_dict[id_index].append(info["read.query_name"])
        df_output = pd.DataFrame(columns=["chrom", "start", "end", "Filter", "Primary_supports", "primary_support",
                                          "Average_Read_Coverage", "mat_gt", "pat_gt", "proband_gt", "total_support_re",
                                          "support_mut", "support_error", "series_support", "failed_support",
                                          "single_support", "shift_err", ])
        for id_index, reads in output_dict.items():
            primary_support = len(reads)
            start, end = id_index.split("_")
            start, end = int(start), int(end)
            this_event_info = {"chrom": chrom, "start": start, "end": end, "Primary_supports": primary_support}
            # print(df_output)
            # filter low support events
            if len(reads) < sample_depth_dict["this_hap"][0] * 0.1:
                this_event_info["Filter"] = "Low_Support_Signal"
                df_output.loc[id_index] = this_event_info
                continue

            # filter abnormal read depth regions
            this_dps = [i.nsegments for i in
                        bam.pileup(contig=chrom,start=start,end=end,min_mapping_quality=mapping_qual_thresh,truncate=True)]
            this_region_depth = 0 if len(this_dps) == 0 else np.nanmean(this_dps)
            this_event_info["Average_Read_Coverage"] = this_region_depth
            if this_region_depth > sample_depth_dict["this_hap"][0] + 3 * sample_depth_dict["this_hap"][1] or \
                    this_region_depth < sample_depth_dict["this_hap"][0] - 3 * sample_depth_dict["this_hap"][1]:
                this_event_info["Filter"] = "Abnormal_Read_Coverage"
                this_event_info["Average_Read_Coverage"] = this_region_depth
                df_output.loc[id_index] = this_event_info
                continue
            # Check  Mendelian inheritance according to bases in reads
            gt_start, gt_end = [], []
            for rec in vcf.fetch(chrom,start - 1,start + 1):
                alts = rec.ref + rec.alts[0]
                gt_start = [alts[i] for i in rec.samples[sample]["GT"]]
            for rec in vcf.fetch(chrom,end - 1,end + 1):
                alts = rec.ref + rec.alts[0]
                gt_end = [alts[i] for i in rec.samples[sample]["GT"]]
            if len(gt_start) != 2 or len(gt_end) != 2:
                this_event_info["Filter"] = "Abnormal_Mutation_Site"
                df_output.loc[id_index] = this_event_info
                continue
            child_hap = [f"{gt_start[0]}_{gt_end[0]}", f"{gt_start[1]}_{gt_end[1]}"]
            expect_rec_hap = [f"{gt_start[0]}_{gt_end[1]}", f"{gt_start[1]}_{gt_end[0]}"]
            for i, j in zip(gt_start,gt_end):
                child_hap.append(f"{i}_{j}")
            proband_gts = {}
            for read in proband_bam.fetch(chrom,start,end):
                if read.reference_start >= start or read.reference_end <= end or read.mapping_quality < mapping_qual_thresh or read.is_secondary or read.is_supplementary: continue
                pairs = {item[1]: item[0] for item in read.get_aligned_pairs(matches_only=True)}
                if start not in pairs or end not in pairs: continue
                start_base = read.query_sequence[pairs[start] - 1]
                end_base = read.query_sequence[pairs[end] - 1]
                gt_read = f"{start_base}_{end_base}"
                if gt_read not in proband_gts: proband_gts[gt_read] = 0
                proband_gts[gt_read] += 1
            pat_gts = {}
            for read in pat_bam.fetch(chrom,start,end):
                if read.reference_start >= start or read.reference_end <= end or read.mapping_quality < mapping_qual_thresh or read.is_secondary or read.is_supplementary: continue
                pairs = {item[1]: item[0] for item in read.get_aligned_pairs(matches_only=True)}
                if start not in pairs or end not in pairs: continue
                start_base = read.query_sequence[pairs[start] - 1]
                end_base = read.query_sequence[pairs[end] - 1]
                gt_read = f"{start_base}_{end_base}"
                if gt_read not in pat_gts: pat_gts[gt_read] = 0
                pat_gts[gt_read] += 1
            mat_gts = {}
            for read in mat_bam.fetch(chrom,start,end):
                if read.reference_start >= start or read.reference_end <= end or read.mapping_quality < mapping_qual_thresh or read.is_secondary or read.is_supplementary: continue
                pairs = {item[1]: item[0] for item in read.get_aligned_pairs(matches_only=True)}
                if start not in pairs or end not in pairs: continue
                start_base = read.query_sequence[pairs[start] - 1]
                end_base = read.query_sequence[pairs[end] - 1]
                gt_read = f"{start_base}_{end_base}"
                if gt_read not in mat_gts: mat_gts[gt_read] = 0
                mat_gts[gt_read] += 1
            pat_gt_list, pat_depth = get_max2_gt(pat_gts,sample_depth_dict["pat"][0])
            mat_gt_list, mat_depth = get_max2_gt(mat_gts,sample_depth_dict["mat"][0])
            proband_gt_list, proband_depth = get_max2_gt(proband_gts,sample_depth_dict["proband"][0])
            this_event_info["mat_gt"] = ",".join(mat_gt_list)
            this_event_info["pat_gt"] = ",".join(pat_gt_list)
            this_event_info["proband_gt"] = ",".join(proband_gt_list)
            if len(proband_gt_list) <= 1 or len(pat_gt_list) <= 1 or len(mat_gt_list) <= 1:
                this_event_info["Filter"] = "Genotype_Unavailable"
                df_output.loc[id_index] = this_event_info
                continue
            else:
                if proband_gt_list[0] in mat_gt_list and proband_gt_list[1] in pat_gt_list:
                    status_mde = False
                elif proband_gt_list[1] in mat_gt_list and proband_gt_list[0] in pat_gt_list:
                    status_mde = False
                else:
                    status_mde = True
            if not status_mde:
                this_event_info["Filter"] = "No_MIE"
                df_output.loc[id_index] = this_event_info
                continue
            support_keep = 0
            support_rec = 0
            support_error = 0
            total_support = 0
            shift_err = 0
            single_support = 0
            series_support = 0
            failed_support = 0
            for read in bam.fetch(chrom,start,end):
                if read.reference_start >= start or read.reference_end <= end or read.mapping_quality < mapping_qual_thresh or read.is_secondary or read.is_supplementary: continue
                pairs = {item[1]: item[0] for item in read.get_aligned_pairs(matches_only=True)}
                if start not in pairs or end not in pairs: continue
                start_base = read.query_sequence[pairs[start] - 1]
                end_base = read.query_sequence[pairs[end] - 1]
                gt_read = f"{start_base}_{end_base}"
                if gt_read in child_hap:
                    support_keep += 1
                elif gt_read in expect_rec_hap:
                    support_rec += 1
                else:
                    support_error += 1
                total_support += 1
                ref_muts_pos = [rec.pos for rec in vcf.fetch(chrom,read.reference_start,read.reference_end)]
                total_mut_num = len(ref_muts_pos)
                vec_reads = []
                vec_hap1 = []
                vec_hap2 = []
                # print(ref_muts_pos)
                mut_pos_in = []
                for pos in ref_muts_pos:
                    if pos not in pairs:
                        # varvec.append("N")
                        continue
                    else:
                        mut_pos_in.append(pos)
                        if (hap1_pots[pos] == 0) and (hap2_pots[pos] == 0):
                            print("error position")
                            print(hap1_pots[pos - 2:pos + 2])
                            print(hap2_pots[pos - 2:pos + 2])
                            exit()
                        read_pos = pairs[pos]
                        read_base = read.query_sequence[read_pos - 1]
                        ref_base = ref.fetch(chrom,pos - 1,pos)
                        if read_base != ref_base:
                            vec_reads.append(1)
                        else:
                            vec_reads.append(0)
                        if hap1_pots[pos] == 1:
                            vec_hap1.append(1)
                            vec_hap2.append(0)
                        else:
                            vec_hap1.append(0)
                            vec_hap2.append(1)
                if len(mut_pos_in) < 2:
                    failed_support += 1
                    continue
                hap1_dis = get_dis(vec_hap1,vec_reads)
                hap2_dis = get_dis(vec_hap2,vec_reads)
                sour_vec = ""
                for i, j in zip(vec_hap1,vec_reads):
                    if i == j:
                        sour_vec += "1"
                    else:
                        sour_vec += "2"
                vec_hap1_str = "".join([f"{i}" for i in vec_hap1])
                vec_hap2_str = "".join([f"{i}" for i in vec_hap2])
                vec_reads_str = "".join([f"{i}" for i in vec_reads])
                frez = sour_vec[0]
                last_i = sour_vec[0]
                shift_index = 0
                for i, index in zip(sour_vec,range(len(sour_vec))):
                    if i == last_i:
                        continue
                    else:
                        shift_index = index
                        frez += i
                        last_i = i
                if len(frez) != 2:
                    shift_err += 1
                    continue
                hap1, hap2 = 0, 0
                for i in sour_vec:
                    if i == "1":
                        hap1 += 1
                    elif i == "2":
                        hap2 += 1
                if hap1 > 1 and hap2 > 1:
                    series_support += 1
                else:
                    single_support += 1
            this_event_info["total_support_re"] = total_support
            this_event_info["support_mut"] = support_rec
            this_event_info["support_error"] = support_error
            this_event_info["series_support"] = series_support
            this_event_info["failed_support"] = failed_support
            this_event_info["single_support"] = single_support  # maybe de novol mutation
            this_event_info["shift_err"] = shift_err  # maybe de novol mutation

            if total_support == 0 or support_rec / (total_support) < 0.8:
                this_event_info["Filter"] = "Maybe_Abnormal_Region"
                print(this_event_info)
                df_output.loc[id_index] = this_event_info
                continue
            if shift_err / total_support > 0.5:
                this_event_info["Filter"] = "Maybe_De_Novol_Mutation"
            else:
                if series_support>shift_err and series_support> sample_depth_dict["this_hap"][0]*0.1:
                    this_event_info["Filter"] = "Meiotic_Recombination"
                else:
                    this_event_info["Filter"] = "Need_more_validation"
            for key, value in this_event_info.items():
                print(key, ": ", value)
            print("+=======================================================")
            df_output.loc[id_index] = this_event_info
        df_output.to_csv(f"{output.read_list}",sep="\t")

rule merge_res:
    input:
        expand(dir_work + "{{family}}/details/{{family}}.{{proband}}.{{hap}}.{{tech}}.phased.byreads.pass.{chrom}.list",chrom=chromsomes)
    output:
        read_list=dir_work + "{family}/{family}.{proband}.{hap}.{tech}.final.merd.tsv",
    run:
        finial_res = pd.DataFrame()
        for item in input:
            this_pd = pd.read_csv(item,index_col=0,sep="\t")
            this_pd = this_pd[this_pd["Filter"] == "Meiotic_Recombination"]
            this_pd.index = [f"{j}_{i}" for i, j in zip(this_pd.index.to_list(),this_pd["chrom"].to_list())]
            finial_res = pd.concat([finial_res, this_pd])
        finial_res.to_csv(f"{output}",sep="\t")


rule tabix:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    run:
        shell("{tabix} {input}")
